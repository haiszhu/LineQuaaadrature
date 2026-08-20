#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#if defined(_WIN32)
  #include <malloc.h>
  #define LQ_ALIGNED_ALLOC(align, size) _aligned_malloc((size), (align))
  #define LQ_ALIGNED_FREE(ptr) _aligned_free(ptr)
#else
  #define LQ_ALIGNED_ALLOC(align, size) aligned_alloc((align), (size))
  #define LQ_ALIGNED_FREE(ptr) free(ptr)
#endif

#if defined(__arm__) || defined(__aarch64__)
  #define LQ_ARCH_ARM
#elif defined(__x86_64__) || defined(_M_X64)
  #if defined(__AVX2__)
    #define LQ_ARCH_AVX2
  #endif
  #if defined(__AVX512F__)
    #define LQ_ARCH_AVX512
  #endif
#endif

#if defined(LQ_ARCH_ARM)
  #include <arm_neon.h>
#elif defined(LQ_ARCH_AVX2) || defined(LQ_ARCH_AVX512)
  #include <immintrin.h>
#endif

#define LQ_RR_MIN 1.0e-300

static void lq_sdlp_scalar(int64_t m, const double *r0, int64_t n,
                           const double *r, const double *rn, const double *w,
                           double *As, double *Ad, int64_t i0)
{
    const double pi4inv = 1.0/(16.0*atan(1.0));
    int64_t i, j;

    for (j = 0; j < n; ++j) {
        double wj = pi4inv*w[j];
        for (i = i0; i < m; ++i) {
            double dx = r0[3*i]   - r[3*j];
            double dy = r0[3*i+1] - r[3*j+1];
            double dz = r0[3*i+2] - r[3*j+2];
            double rr = dx*dx + dy*dy + dz*dz;
            if (rr > LQ_RR_MIN) {
                double sj = wj*(1.0/sqrt(rr));
                As[i+j*m] = sj;
                Ad[i+j*m] = sj*(dx*rn[3*j]+dy*rn[3*j+1]+dz*rn[3*j+2])/rr;
            } else { As[i+j*m] = 0.0;  Ad[i+j*m] = 0.0; }
        }
    }
}

static double *lq_soa_alloc(int64_t mv, size_t align, const double *r0,
                            double **r0x, double **r0y, double **r0z)
{
    double *buf;
    int64_t i;
    if (mv <= 0) { *r0x = *r0y = *r0z = NULL;  return NULL; }
    buf = (double*)LQ_ALIGNED_ALLOC(align, 3*mv*sizeof(double));
    *r0x = buf;  *r0y = buf + mv;  *r0z = buf + 2*mv;
    for (i = 0; i < mv; ++i) {
        (*r0x)[i] = r0[3*i];  (*r0y)[i] = r0[3*i+1];  (*r0z)[i] = r0[3*i+2];
    }
    return buf;
}

void lq_csimd128lap3dsdlpmat_c_(int64_t *M, const double *r0, int64_t *N,
                                const double *r, const double *rn,
                                const double *w, double *As, double *Ad)
{
#ifdef LQ_ARCH_ARM
    const int64_t m = *M, n = *N;
    const double  pi4inv = 1.0/(16.0*atan(1.0));
    const int64_t m2 = (m/2)*2;
    double *r0x, *r0y, *r0z;
    double *buf = lq_soa_alloc(m2, 16, r0, &r0x, &r0y, &r0z);
    const float64x2_t one  = vdupq_n_f64(1.0);
    const float64x2_t zero = vdupq_n_f64(0.0);
    const float64x2_t rmin = vdupq_n_f64(LQ_RR_MIN);
    int64_t i, j;

    for (j = 0; j < n; ++j) {
        const float64x2_t rx_v  = vdupq_n_f64(r[3*j]);
        const float64x2_t ry_v  = vdupq_n_f64(r[3*j+1]);
        const float64x2_t rz_v  = vdupq_n_f64(r[3*j+2]);
        const float64x2_t rnx_v = vdupq_n_f64(rn[3*j]);
        const float64x2_t rny_v = vdupq_n_f64(rn[3*j+1]);
        const float64x2_t rnz_v = vdupq_n_f64(rn[3*j+2]);
        const float64x2_t w_v   = vdupq_n_f64(pi4inv*w[j]);
        double *Asj = As + j*m, *Adj = Ad + j*m;

        for (i = 0; i <= m2 - 2; i += 2) {
            float64x2_t dx_v = vsubq_f64(vld1q_f64(&r0x[i]), rx_v);
            float64x2_t dy_v = vsubq_f64(vld1q_f64(&r0y[i]), ry_v);
            float64x2_t dz_v = vsubq_f64(vld1q_f64(&r0z[i]), rz_v);
            float64x2_t rr_v = vaddq_f64(vaddq_f64(vmulq_f64(dx_v, dx_v),
                                                   vmulq_f64(dy_v, dy_v)),
                                         vmulq_f64(dz_v, dz_v));
            float64x2_t rdotrn_v = vaddq_f64(vaddq_f64(vmulq_f64(dx_v, rnx_v),
                                                       vmulq_f64(dy_v, rny_v)),
                                             vmulq_f64(dz_v, rnz_v));
            uint64x2_t  ok_v = vcgtq_f64(rr_v, rmin);
            float64x2_t rs_v = vbslq_f64(ok_v, rr_v, one);
            float64x2_t ri_v = vdivq_f64(one, vsqrtq_f64(rs_v));
            float64x2_t slpij_v = vbslq_f64(ok_v, vmulq_f64(w_v, ri_v), zero);
            float64x2_t dlpij_v = vdivq_f64(vmulq_f64(slpij_v, rdotrn_v), rs_v);
            vst1q_f64(&Asj[i], slpij_v);
            vst1q_f64(&Adj[i], dlpij_v);
        }
    }
    lq_sdlp_scalar(m, r0, n, r, rn, w, As, Ad, m2);
    LQ_ALIGNED_FREE(buf);
#else
    lq_sdlp_scalar(*M, r0, *N, r, rn, w, As, Ad, 0);
#endif
}

void lq_csimd256lap3dsdlpmat_c_(int64_t *M, const double *r0, int64_t *N,
                                const double *r, const double *rn,
                                const double *w, double *As, double *Ad)
{
#ifdef LQ_ARCH_AVX2
    const int64_t m = *M, n = *N;
    const double  pi4inv = 1.0/(16.0*atan(1.0));
    const int64_t m4 = (m/4)*4;
    double *r0x, *r0y, *r0z;
    double *buf = lq_soa_alloc(m4, 32, r0, &r0x, &r0y, &r0z);
    const __m256d half  = _mm256_set1_pd(0.5);
    const __m256d three = _mm256_set1_pd(3.0);
    const __m256d zero  = _mm256_setzero_pd();
    const __m256d rmin  = _mm256_set1_pd(LQ_RR_MIN);
    int64_t i, j;

    for (j = 0; j < n; ++j) {
        const __m256d rx_v  = _mm256_set1_pd(r[3*j]);
        const __m256d ry_v  = _mm256_set1_pd(r[3*j+1]);
        const __m256d rz_v  = _mm256_set1_pd(r[3*j+2]);
        const __m256d rnx_v = _mm256_set1_pd(rn[3*j]);
        const __m256d rny_v = _mm256_set1_pd(rn[3*j+1]);
        const __m256d rnz_v = _mm256_set1_pd(rn[3*j+2]);
        const __m256d w_v   = _mm256_set1_pd(pi4inv*w[j]);
        double *Asj = As + j*m, *Adj = Ad + j*m;

        for (i = 0; i <= m4 - 4; i += 4) {
            __m256d dx_v = _mm256_sub_pd(_mm256_load_pd(&r0x[i]), rx_v);
            __m256d dy_v = _mm256_sub_pd(_mm256_load_pd(&r0y[i]), ry_v);
            __m256d dz_v = _mm256_sub_pd(_mm256_load_pd(&r0z[i]), rz_v);
            __m256d rr_v = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx_v, dx_v),
                                                       _mm256_mul_pd(dy_v, dy_v)),
                                         _mm256_mul_pd(dz_v, dz_v));
            __m256d rdotrn_v = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx_v, rnx_v),
                                                           _mm256_mul_pd(dy_v, rny_v)),
                                             _mm256_mul_pd(dz_v, rnz_v));
            __m256d ok_v = _mm256_cmp_pd(rr_v, rmin, _CMP_GT_OQ);
            __m256d rs_v = _mm256_blendv_pd(_mm256_set1_pd(1.0), rr_v, ok_v);
            __m256d ri_v = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(rs_v)));
            __m256d muls = _mm256_mul_pd(_mm256_mul_pd(rs_v, ri_v), ri_v);
            ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls));
            muls = _mm256_mul_pd(_mm256_mul_pd(rs_v, ri_v), ri_v);
            ri_v = _mm256_mul_pd(_mm256_mul_pd(half, ri_v), _mm256_sub_pd(three, muls));
            __m256d slpij_v = _mm256_blendv_pd(zero, _mm256_mul_pd(w_v, ri_v), ok_v);
            __m256d ri3_v   = _mm256_mul_pd(_mm256_mul_pd(ri_v, ri_v), ri_v);
            __m256d dlpij_v = _mm256_blendv_pd(zero,
                              _mm256_mul_pd(_mm256_mul_pd(w_v, rdotrn_v), ri3_v), ok_v);
            _mm256_storeu_pd(&Asj[i], slpij_v);
            _mm256_storeu_pd(&Adj[i], dlpij_v);
        }
    }
    lq_sdlp_scalar(m, r0, n, r, rn, w, As, Ad, m4);
    LQ_ALIGNED_FREE(buf);
#else
    lq_sdlp_scalar(*M, r0, *N, r, rn, w, As, Ad, 0);
#endif
}

void lq_csimd512lap3dsdlpmat_c_(int64_t *M, const double *r0, int64_t *N,
                                const double *r, const double *rn,
                                const double *w, double *As, double *Ad)
{
#ifdef LQ_ARCH_AVX512
    const int64_t m = *M, n = *N;
    const double  pi4inv = 1.0/(16.0*atan(1.0));
    const int64_t m8 = (m/8)*8;
    double *r0x, *r0y, *r0z;
    double *buf = lq_soa_alloc(m8, 64, r0, &r0x, &r0y, &r0z);
    const __m512d half  = _mm512_set1_pd(0.5);
    const __m512d three = _mm512_set1_pd(3.0);
    const __m512d one   = _mm512_set1_pd(1.0);
    const __m512d zero  = _mm512_setzero_pd();
    const __m512d rmin  = _mm512_set1_pd(LQ_RR_MIN);
    int64_t i, j;

    for (j = 0; j < n; ++j) {
        const __m512d rx_v  = _mm512_set1_pd(r[3*j]);
        const __m512d ry_v  = _mm512_set1_pd(r[3*j+1]);
        const __m512d rz_v  = _mm512_set1_pd(r[3*j+2]);
        const __m512d rnx_v = _mm512_set1_pd(rn[3*j]);
        const __m512d rny_v = _mm512_set1_pd(rn[3*j+1]);
        const __m512d rnz_v = _mm512_set1_pd(rn[3*j+2]);
        const __m512d w_v   = _mm512_set1_pd(pi4inv*w[j]);
        double *Asj = As + j*m, *Adj = Ad + j*m;

        for (i = 0; i <= m8 - 8; i += 8) {
            __m512d dx_v = _mm512_sub_pd(_mm512_load_pd(&r0x[i]), rx_v);
            __m512d dy_v = _mm512_sub_pd(_mm512_load_pd(&r0y[i]), ry_v);
            __m512d dz_v = _mm512_sub_pd(_mm512_load_pd(&r0z[i]), rz_v);
            __m512d rr_v = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(dx_v, dx_v),
                                                       _mm512_mul_pd(dy_v, dy_v)),
                                         _mm512_mul_pd(dz_v, dz_v));
            __m512d rdotrn_v = _mm512_add_pd(_mm512_add_pd(_mm512_mul_pd(dx_v, rnx_v),
                                                           _mm512_mul_pd(dy_v, rny_v)),
                                             _mm512_mul_pd(dz_v, rnz_v));
            __mmask8 ok_v = _mm512_cmp_pd_mask(rr_v, rmin, _CMP_GT_OQ);
            __m512d  rs_v = _mm512_mask_blend_pd(ok_v, one, rr_v);
            __m512d  ri_v = _mm512_cvtps_pd(_mm256_rsqrt_ps(_mm512_cvtpd_ps(rs_v)));
            __m512d  muls = _mm512_mul_pd(_mm512_mul_pd(rs_v, ri_v), ri_v);
            ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls));
            muls = _mm512_mul_pd(_mm512_mul_pd(rs_v, ri_v), ri_v);
            ri_v = _mm512_mul_pd(_mm512_mul_pd(half, ri_v), _mm512_sub_pd(three, muls));
            __m512d slpij_v = _mm512_mask_blend_pd(ok_v, zero, _mm512_mul_pd(w_v, ri_v));
            __m512d ri3_v   = _mm512_mul_pd(_mm512_mul_pd(ri_v, ri_v), ri_v);
            __m512d dlpij_v = _mm512_mask_blend_pd(ok_v, zero,
                              _mm512_mul_pd(_mm512_mul_pd(w_v, rdotrn_v), ri3_v));
            _mm512_storeu_pd(&Asj[i], slpij_v);
            _mm512_storeu_pd(&Adj[i], dlpij_v);
        }
    }
    lq_sdlp_scalar(m, r0, n, r, rn, w, As, Ad, m8);
    LQ_ALIGNED_FREE(buf);
#else
    lq_sdlp_scalar(*M, r0, *N, r, rn, w, As, Ad, 0);
#endif
}
