package ckks_fv

import (
	"fmt"
	"math"

	"github.com/ldsec/lattigo/v2/ckks/bettersine"
	"github.com/ldsec/lattigo/v2/utils"
	"github.com/ldsec/lattigo/v2/ring"

	"reflect"
	"unsafe"
)

// ShallowCopy 自分で追加
/*func (hbtp *HalfBootstrapper) ShallowCopy() *HalfBootstrapper {
	evalCopy := hbtp.ckksEvaluator.ShallowCopy().(*ckksEvaluator)
	evalCopy.ctxpool = NewCiphertextCKKS(hbtp.params, 1, hbtp.params.MaxLevel(), 0)

	var pDFTInvCopy []*PtDiagMatrix
    if hbtp.pDFTInvWithoutRepack != nil {
        pDFTInvCopy = make([]*PtDiagMatrix, len(hbtp.pDFTInvWithoutRepack))
        for i, mat := range hbtp.pDFTInvWithoutRepack {
			if mat != nil {
				clonedMat := &PtDiagMatrix{
					LogSlots:   mat.LogSlots,
					N1:         mat.N1,
					Level:      mat.Level,
					Scale:      mat.Scale,
					Vec:        make(map[int][2]*ring.Poly, len(mat.Vec)), // 新しいマップを確保
					naive:      mat.naive,
					isGaussian: mat.isGaussian,
				}
				for k, v := range mat.Vec {
					var poly0, poly1 *ring.Poly
					if v[0] != nil {
						poly0 = v[0].CopyNew()
					}
					if v[1] != nil {
						poly1 = v[1].CopyNew()
					}
					clonedMat.Vec[k] = [2]*ring.Poly{poly0, poly1}
				}
				pDFTInvCopy[i] = clonedMat
            }
        }
    }

	return &HalfBootstrapper{
		ckksEvaluator:      evalCopy,
		HalfBootParameters: hbtp.HalfBootParameters,
		BootstrappingKey: &BootstrappingKey{hbtp.BootstrappingKey.Rlk, hbtp.BootstrappingKey.Rtks},

		params: hbtp.params.Copy(),
		dslots: hbtp.dslots,
		logdslots: hbtp.logdslots,

		encoder: NewCKKSEncoder(hbtp.params),

		prescale: hbtp.prescale,
		postscale: hbtp.postscale,
		sinescale: hbtp.sinescale,
		sqrt2pi: hbtp.sqrt2pi,
		scFac: hbtp.scFac,
		sineEvalPoly: hbtp.sineEvalPoly,
		arcSinePoly: hbtp.arcSinePoly,

		coeffsToSlotsDiffScale: hbtp.coeffsToSlotsDiffScale,
		diffScaleAfterSineEval: hbtp.diffScaleAfterSineEval,
		pDFTInvWithoutRepack: pDFTInvCopy,

		rotKeyIndex: hbtp.rotKeyIndex,
	}
}*/

// ShallowCopy 自分で追加（完全監査クリア版）
func (hbtp *HalfBootstrapper) ShallowCopy() *HalfBootstrapper {
	var err error

	if hbtp == nil {
		return nil
	}

	baseCopy := *hbtp.ckksEvaluator.ckksEvaluatorBase

	// ckksEvaluatorBaseのringQ, ringPを新調
	if hbtp.ckksEvaluatorBase.ringQ != nil {
		baseCopy.ringQ, err = ring.NewRing(hbtp.ckksEvaluatorBase.ringQ.N, hbtp.ckksEvaluatorBase.ringQ.Modulus)
		if err != nil {
			panic(fmt.Errorf("ShallowCopy: failed to recreate ringQ: %w", err))
		}
	}
	if hbtp.ckksEvaluatorBase.ringP != nil {
		baseCopy.ringP, err = ring.NewRing(hbtp.ckksEvaluatorBase.ringP.N, hbtp.ckksEvaluatorBase.ringP.Modulus)
		if err != nil {
			panic(fmt.Errorf("ShallowCopy: failed to recreate ringP: %w", err))
		}
	}

	// baseconverterも新調
	var baseconverterCopy *ring.FastBasisExtender
	if hbtp.params.PiCount() != 0 {
		baseconverterCopy = ring.NewFastBasisExtender(baseCopy.ringQ, baseCopy.ringP)
	}

	buffCopy := newCKKSEvaluatorBuffers(&baseCopy)
	buffCopy.ctxpool = NewCiphertextCKKS(hbtp.params, 1, hbtp.params.MaxLevel(), 0)

	evalCopy := &ckksEvaluator{
		ckksEvaluatorBase:    &baseCopy,
		ckksEvaluatorBuffers: buffCopy,
		rlk:                  hbtp.ckksEvaluator.rlk,
		rtks:                 hbtp.ckksEvaluator.rtks,
		permuteNTTIndex:      hbtp.ckksEvaluator.permuteNTTIndex,
		baseconverter:        baseconverterCopy,
	}

	return &HalfBootstrapper{
		ckksEvaluator:      evalCopy,
		HalfBootParameters: hbtp.HalfBootParameters,
		//BootstrappingKey:   &BootstrappingKey{Rlk: hbtp.BootstrappingKey.Rlk, Rtks: hbtp.BootstrappingKey.Rtks},
		BootstrappingKey: 	hbtp.BootstrappingKey,

		params:    hbtp.params,
		dslots:    hbtp.dslots,
		logdslots: hbtp.logdslots,

		encoder: NewCKKSEncoder(hbtp.params),

		prescale:     hbtp.prescale,
		postscale:    hbtp.postscale,
		sinescale:    hbtp.sinescale,
		sqrt2pi:      hbtp.sqrt2pi,
		scFac:        hbtp.scFac,
		sineEvalPoly: hbtp.sineEvalPoly,
		arcSinePoly:  hbtp.arcSinePoly,

		coeffsToSlotsDiffScale: hbtp.coeffsToSlotsDiffScale,
		diffScaleAfterSineEval: hbtp.diffScaleAfterSineEval,
		//pDFTInvWithoutRepack:   pDFTInvCopy,
		pDFTInvWithoutRepack:   hbtp.pDFTInvWithoutRepack,

		rotKeyIndex: hbtp.rotKeyIndex,
	}
}



// HalfBootstrapper is a struct to stores a memory pool the plaintext matrices
// the polynomial approximation and the keys for the half-bootstrapping.
type HalfBootstrapper struct {
	*ckksEvaluator
	HalfBootParameters
	*BootstrappingKey
	params *Parameters

	dslots    int // Number of plaintext slots after the re-encoding
	logdslots int

	encoder CKKSEncoder // Encoder

	prescale     float64                 // Q[0]/(Q[0]/|m|)
	postscale    float64                 // Qi sineeval/(Q[0]/|m|)
	sinescale    float64                 // Qi sineeval
	sqrt2pi      float64                 // (1/2pi)^{-2^r}
	scFac        float64                 // 2^{r}
	sineEvalPoly *ChebyshevInterpolation // Coefficients of the Chebyshev Interpolation of sin(2*pi*x) or cos(2*pi*x/r)
	arcSinePoly  *Poly                   // Coefficients of the Taylor series of arcsine(x)

	coeffsToSlotsDiffScale complex128      // Matrice rescaling
	diffScaleAfterSineEval float64         // Matrice rescaling
	pDFTInvWithoutRepack   []*PtDiagMatrix // Matrice vectors

	rotKeyIndex []int // a list of the required rotation keys
}

// NewHalfBootstrapper creates a new HalfBootstrapper.
func NewHalfBootstrapper(params *Parameters, hbtpParams *HalfBootParameters, btpKey BootstrappingKey) (hbtp *HalfBootstrapper, err error) {

	if hbtpParams.SinType == SinType(Sin) && hbtpParams.SinRescal != 0 {
		return nil, fmt.Errorf("cannot use double angle formul for SinType = Sin -> must use SinType = Cos")
	}

	hbtp = newHalfBootstrapper(params, hbtpParams)

	hbtp.BootstrappingKey = &BootstrappingKey{btpKey.Rlk, btpKey.Rtks}
	if err = hbtp.CheckKeys(); err != nil {
		return nil, fmt.Errorf("invalid bootstrapping key: %w", err)
	}
	hbtp.ckksEvaluator = hbtp.ckksEvaluator.WithKey(EvaluationKey{btpKey.Rlk, btpKey.Rtks}).(*ckksEvaluator)

	return hbtp, nil
}

// newHalfBootstrapper is a constructor of "dummy" half-bootstrapper to enable the generation of bootstrapping-related constants
// without providing a bootstrapping key. To be replaced by a propper factorization of the bootstrapping pre-computations.
func newHalfBootstrapper(params *Parameters, hbtpParams *HalfBootParameters) (hbtp *HalfBootstrapper) {
	hbtp = new(HalfBootstrapper)

	hbtp.params = params.Copy()
	hbtp.HalfBootParameters = *hbtpParams.Copy()

	hbtp.dslots = params.Slots()
	hbtp.logdslots = params.LogSlots()
	if params.logSlots < params.MaxLogSlots() {
		hbtp.dslots <<= 1
		hbtp.logdslots++
	}

	hbtp.prescale = math.Exp2(math.Round(math.Log2(float64(params.qi[0]) / hbtp.MessageRatio)))
	hbtp.sinescale = math.Exp2(math.Round(math.Log2(hbtp.SineEvalModuli.ScalingFactor)))
	hbtp.postscale = hbtp.sinescale / hbtp.MessageRatio

	hbtp.encoder = NewCKKSEncoder(params)
	hbtp.ckksEvaluator = NewCKKSEvaluator(params, EvaluationKey{}).(*ckksEvaluator) // creates an evaluator without keys for genDFTMatrices

	hbtp.genSinePoly()
	hbtp.genDFTMatrices()

	hbtp.ctxpool = NewCiphertextCKKS(params, 1, params.MaxLevel(), 0)

	return hbtp
}

// CheckKeys checks if all the necessary keys are present
func (hbtp *HalfBootstrapper) CheckKeys() (err error) {

	if hbtp.Rlk == nil {
		return fmt.Errorf("relinearization key is nil")
	}

	if hbtp.Rtks == nil {
		return fmt.Errorf("rotation key is nil")
	}

	rotMissing := []int{}
	for _, i := range hbtp.rotKeyIndex {
		galEl := hbtp.params.GaloisElementForColumnRotationBy(int(i))
		if _, generated := hbtp.Rtks.Keys[galEl]; !generated {
			rotMissing = append(rotMissing, i)
		}
	}

	if len(rotMissing) != 0 {
		return fmt.Errorf("rotation key(s) missing: %d", rotMissing)
	}

	return nil
}

func (hbtp *HalfBootstrapper) genDFTMatrices() {

	a := real(hbtp.sineEvalPoly.a)
	b := real(hbtp.sineEvalPoly.b)
	n := float64(hbtp.params.N())
	qDiff := float64(hbtp.params.qi[0]) / math.Exp2(math.Round(math.Log2(float64(hbtp.params.qi[0]))))

	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + evantual scaling factor for the double angle formula
	hbtp.coeffsToSlotsDiffScale = complex(math.Pow(2.0/((b-a)*n*hbtp.scFac*qDiff), 1.0/float64(hbtp.CtSDepth(false))), 0)

	// Rescaling factor to set the final ciphertext to the desired scale
	hbtp.diffScaleAfterSineEval = (qDiff * hbtp.params.scale) / hbtp.postscale

	// CoeffsToSlotsWithoutRepack vectors
	hbtp.pDFTInvWithoutRepack = hbtp.HalfBootParameters.GenCoeffsToSlotsMatrixWithoutRepack(hbtp.coeffsToSlotsDiffScale, hbtp.encoder)

	// List of the rotation key values to needed for the bootstrapp
	hbtp.rotKeyIndex = []int{}

	//SubSum rotation needed X -> Y^slots rotations
	for i := hbtp.params.logSlots; i < hbtp.params.MaxLogSlots(); i++ {
		if !utils.IsInSliceInt(1<<i, hbtp.rotKeyIndex) {
			hbtp.rotKeyIndex = append(hbtp.rotKeyIndex, 1<<i)
		}
	}

	// Coeffs to Slots rotations
	for _, pVec := range hbtp.pDFTInvWithoutRepack {
		hbtp.rotKeyIndex = AddMatrixRotToList(pVec, hbtp.rotKeyIndex, hbtp.params.Slots(), false)
	}
}

func (hbtp *HalfBootstrapper) genSinePoly() {

	K := int(hbtp.SinRange)
	deg := int(hbtp.SinDeg)
	hbtp.scFac = float64(int(1 << hbtp.SinRescal))

	if hbtp.ArcSineDeg > 0 {
		hbtp.sqrt2pi = 1.0

		coeffs := make([]complex128, hbtp.ArcSineDeg+1)

		coeffs[1] = 0.15915494309189535

		for i := 3; i < hbtp.ArcSineDeg+1; i += 2 {

			coeffs[i] = coeffs[i-2] * complex(float64(i*i-4*i+4)/float64(i*i-i), 0)

		}

		hbtp.arcSinePoly = NewPoly(coeffs)

	} else {
		hbtp.sqrt2pi = math.Pow(0.15915494309189535, 1.0/hbtp.scFac)
	}

	if hbtp.SinType == Sin {

		hbtp.sineEvalPoly = Approximate(sin2pi2pi, -complex(float64(K)/hbtp.scFac, 0), complex(float64(K)/hbtp.scFac, 0), deg)

	} else if hbtp.SinType == Cos1 {

		hbtp.sineEvalPoly = new(ChebyshevInterpolation)

		hbtp.sineEvalPoly.coeffs = bettersine.Approximate(K, deg, hbtp.MessageRatio, int(hbtp.SinRescal))

		hbtp.sineEvalPoly.maxDeg = hbtp.sineEvalPoly.Degree()
		hbtp.sineEvalPoly.a = complex(float64(-K)/hbtp.scFac, 0)
		hbtp.sineEvalPoly.b = complex(float64(K)/hbtp.scFac, 0)
		hbtp.sineEvalPoly.lead = true

	} else if hbtp.SinType == Cos2 {

		hbtp.sineEvalPoly = Approximate(cos2pi, -complex(float64(K)/hbtp.scFac, 0), complex(float64(K)/hbtp.scFac, 0), deg)

	} else {
		panic("Bootstrapper -> invalid sineType")
	}

	for i := range hbtp.sineEvalPoly.coeffs {
		hbtp.sineEvalPoly.coeffs[i] *= complex(hbtp.sqrt2pi, 0)
	}
}







// HbtpCrossSpy は、構造体の最深部の「値」を記録するスパイ
type HbtpCrossSpy struct {
	snapshot map[string]interface{}
}

// NewHbtpCrossSpy は、HalfBootstrapper の全最深部データを記録する
func NewHbtpCrossSpy(hbtp *HalfBootstrapper) *HbtpCrossSpy {
	spy := &HbtpCrossSpy{snapshot: make(map[string]interface{})}
	if hbtp == nil {
		return spy
	}
	spy.scan(reflect.ValueOf(hbtp).Elem(), "HalfBootstrapper")
	return spy
}

func (s *HbtpCrossSpy) scan(v reflect.Value, path string) {
	if !v.IsValid() { return }

	for v.Kind() == reflect.Ptr || v.Kind() == reflect.Interface {
		if v.IsNil() {
			s.snapshot[path+" (nil)"] = nil
			return
		}
		v = v.Elem()
	}

	switch v.Kind() {
	case reflect.Struct:
		t := v.Type()
		for i := 0; i < v.NumField(); i++ {
			f := v.Field(i)
			fName := t.Field(i).Name
			// 小文字フィールドも強制展開
			f = reflect.NewAt(f.Type(), unsafe.Pointer(f.UnsafeAddr())).Elem()
			s.scan(f, path+"."+fName)
		}
	case reflect.Slice, reflect.Array:
		s.snapshot[path+"_len"] = v.Len()
		maxLen := v.Len()
		if maxLen > 8 { maxLen = 8 }
		for i := 0; i < maxLen; i++ {
			s.scan(v.Index(i), fmt.Sprintf("%s[%d]", path, i))
		}
	case reflect.Map:
		s.snapshot[path+"_len"] = v.Len()
		iter := v.MapRange()
		count := 0
		for iter.Next() {
			if count > 5 { break }
			s.scan(iter.Value(), fmt.Sprintf("%s[map_key_%v]", path, iter.Key().Interface()))
			count++
		}
	default:
		if v.CanInterface() {
			s.snapshot[path] = v.Interface()
		}
	}
}

// CheckContamination は、1周目の計算のせいで、2周目（触っていないはずの側）に変化が飛び火したかを暴く
func (s *HbtpCrossSpy) CheckContamination(afterClean *HbtpCrossSpy) {
	fmt.Println("\n🚨 ===== [CROSS-INSTANCE CONTAMINATION REPORT] =====")
	contaminated := false

	for path, bVal := range s.snapshot {
		aVal, exists := afterClean.snapshot[path]
		if !exists { continue }

		if !reflect.DeepEqual(bVal, aVal) {
			fmt.Printf("❌ [POINTER LEAK DETECTED] フィールド '%s' が裏で繋がっています！\n", path)
			fmt.Printf("   ├─ 計算前の状態: %v\n", bVal)
			fmt.Printf("   └─ 1周目計算後の状態: %v\n", aVal)
			contaminated = true
		}
	}

	if !contaminated {
		fmt.Println("✅ [PERFECT ISOLATION] ShallowCopy されたインスタンス同士は、裏で一切繋がっていません！完全隔離されています。")
	}
	fmt.Println("=====================================================\n")
}