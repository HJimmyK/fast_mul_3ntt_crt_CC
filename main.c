// MIT License
// 
// Copyright (c) 2025 Jecricho Knox
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <time.h>
#include "core.h"
#include <stdlib.h>
#include <stdio.h>

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

#define ALIGNED_MALLOC(_ptr, _type, _len)                                                             \
    do {                                                                                              \
        size_t elem_size__ = sizeof(_type);                                                           \
        size_t total_size__ = (_len) * elem_size__;                                                   \
        (_ptr) = (total_size__ == 0) ? NULL : (_type*)_aligned_malloc(total_size__, CACHE_LINE_SIZE); \
    } while (0)

#define ALIGNED_FREE(_ptr)       \
    do {                         \
        if ((_ptr) != NULL) {    \
            _aligned_free(_ptr); \
            (_ptr) = NULL;       \
        }                        \
    } while (0)

#define define_convolution_rec(_i)                                                                                    \
    void convolution_rec_##_i(mont64* in1, mont64* in2, mont64* out, ntt_short* table, size_t ntt_len, bool norm) {   \
        assert(in1 != NULL && in2 != NULL && out != NULL && table != NULL);                                           \
        assert(in1 != in2);                                                                                           \
        if (ntt_len <= long_threshold) {                                                                              \
            dif_func(in1, table, ntt_len, _i);                                                                        \
            dif_func(in2, table, ntt_len, _i);                                                                        \
            if (norm) {                                                                                               \
                mont64 inv_len = ntt_len;                                                                             \
                _mont_tomont_func(inv_len, _i);                                                                       \
                inv_len = _mont_qpow_func_name(_i)(inv_len, ((g_mod(_i)) - 2));                                       \
                for (size_t ii = 0; ii < ntt_len; ii++) {                                                             \
                    _mont_mul_func(out[ii], in1[ii], in2[ii], _i);                                                    \
                    _mont_mulinto_func(out[ii], inv_len, _i);                                                         \
                }                                                                                                     \
            } else {                                                                                                  \
                for (size_t ii = 0; ii < ntt_len; ii++) {                                                             \
                    _mont_mul_func(out[ii], in1[ii], in2[ii], _i);                                                    \
                }                                                                                                     \
            }                                                                                                         \
            idit_func(out, table, ntt_len, _i);                                                                       \
            return;                                                                                                   \
        }                                                                                                             \
        const size_t quarter_len = ntt_len / 4;                                                                       \
        mont64 unit_omega1 = g_root(_i);                                                                              \
        _mont_tomont_func(unit_omega1, _i);                                                                           \
        unit_omega1 = _mont_qpow_func_name(_i)(unit_omega1, (g_mod(_i) - 1) / ntt_len);                               \
        mont64 unit_omega3 = _mont_qpow_func_name(_i)(unit_omega1, 3);                                                \
        mont64 omega1 = g_one(_i), omega3 = g_one(_i);                                                                \
        for (size_t ii = 0; ii < quarter_len; ii++) {                                                                 \
            mont64 temp0 = in1[ii], temp1 = in1[quarter_len + ii];                                                    \
            mont64 temp2 = in1[quarter_len * 2 + ii], temp3 = in1[quarter_len * 3 + ii];                              \
            _dif_butterfly244(temp0, temp1, temp2, temp3, _i);                                                        \
            in1[ii] = temp0, in1[quarter_len + ii] = temp1;                                                           \
            _mont_mul_func(in1[quarter_len * 2 + ii], temp2, omega1, _i);                                             \
            _mont_mul_func(in1[quarter_len * 3 + ii], temp3, omega3, _i);                                             \
            temp0 = in2[ii], temp1 = in2[quarter_len + ii];                                                           \
            temp2 = in2[quarter_len * 2 + ii], temp3 = in2[quarter_len * 3 + ii];                                     \
            _dif_butterfly244(temp0, temp1, temp2, temp3, _i);                                                        \
            in2[ii] = temp0, in2[quarter_len + ii] = temp1;                                                           \
            _mont_mul_func(in2[quarter_len * 2 + ii], temp2, omega1, _i);                                             \
            _mont_mul_func(in2[quarter_len * 3 + ii], temp3, omega3, _i);                                             \
            _mont_mulinto_func(omega1, unit_omega1, _i);                                                              \
            _mont_mulinto_func(omega3, unit_omega3, _i);                                                              \
        }                                                                                                             \
        convolution_rec_##_i(in1, in2, out, table, ntt_len / 2, false);                                               \
        convolution_rec_##_i(in1 + quarter_len * 2, in2 + quarter_len * 2, out + quarter_len * 2, table, ntt_len / 4, \
                             false);                                                                                  \
        convolution_rec_##_i(in1 + quarter_len * 3, in2 + quarter_len * 3, out + quarter_len * 3, table, ntt_len / 4, \
                             false);                                                                                  \
        unit_omega1 = g_rootinv(_i);                                                                                  \
        _mont_tomont_func(unit_omega1, _i);                                                                           \
        unit_omega1 = _mont_qpow_func_name(_i)(unit_omega1, (g_mod(_i) - 1) / ntt_len);                               \
        unit_omega3 = _mont_qpow_func_name(_i)(unit_omega1, 3);                                                       \
        if (norm) {                                                                                                   \
            mont64 inv_len = ntt_len;                                                                                 \
            _mont_tomont_func(inv_len, _i);                                                                           \
            inv_len = _mont_qpow_func_name(_i)(inv_len, (g_mod(_i) - 2));                                             \
            omega1 = inv_len, omega3 = inv_len;                                                                       \
            for (size_t ii = 0; ii < quarter_len; ii++) {                                                             \
                mont64 temp0, temp1, temp2, temp3;                                                                    \
                _mont_mul_func(temp0, out[ii], inv_len, _i);                                                          \
                _mont_mul_func(temp1, out[quarter_len + ii], inv_len, _i);                                            \
                _mont_mul_func(temp2, out[quarter_len * 2 + ii], omega1, _i);                                         \
                _mont_mul_func(temp3, out[quarter_len * 3 + ii], omega3, _i);                                         \
                _idit_butterfly244(temp0, temp1, temp2, temp3, _i);                                                   \
                out[ii] = temp0, out[quarter_len + ii] = temp1;                                                       \
                out[quarter_len * 2 + ii] = temp2, out[quarter_len * 3 + ii] = temp3;                                 \
                _mont_mulinto_func(omega1, unit_omega1, _i);                                                          \
                _mont_mulinto_func(omega3, unit_omega3, _i);                                                          \
            }                                                                                                         \
        } else {                                                                                                      \
            omega1 = g_one(_i), omega3 = g_one(_i);                                                                   \
            for (size_t ii = 0; ii < quarter_len; ii++) {                                                             \
                mont64 temp0 = out[ii], temp1 = out[quarter_len + ii], temp2, temp3;                                  \
                _mont_mul_func(temp2, out[quarter_len * 2 + ii], omega1, _i);                                         \
                _mont_mul_func(temp3, out[quarter_len * 3 + ii], omega3, _i);                                         \
                _idit_butterfly244(temp0, temp1, temp2, temp3, _i);                                                   \
                out[ii] = temp0, out[quarter_len + ii] = temp1;                                                       \
                out[quarter_len * 2 + ii] = temp2, out[quarter_len * 3 + ii] = temp3;                                 \
                _mont_mulinto_func(omega1, unit_omega1, _i);                                                          \
                _mont_mulinto_func(omega3, unit_omega3, _i);                                                          \
            }                                                                                                         \
        }                                                                                                             \
    }

define_convolution_rec(1) define_convolution_rec(2) define_convolution_rec(3)

#define convolution_rec_func(in1, in2, out, table, ntt_len, _i) convolution_rec_##_i(in1, in2, out, table, ntt_len, true)

u64 int_ceil2(u64 n) {
    const int bits = 64;
    n--;
    for (int i = 1; i < bits; i *= 2) {
        n |= (n >> i);
    }
    return n + 1;
}

//static mont64 g_a[long_threshold];
//static mont64 g_b[long_threshold];
double time_ntt[251][5];

void abs_mul64(u64* in1, u64 len1, u64* in2, u64 len2, u64* out, uint64_t index) {
    assert(in1 != NULL && in2 != NULL && out != NULL);
    u64 out_len = len1 + len2, conv_len = out_len - 1;
    u64 ntt_len = int_ceil2(conv_len);

    mont64* buf1_mont, *buf2_mont, *buf3_mont, *tmp_mont;
    ALIGNED_MALLOC(buf1_mont, mont64, ntt_len);
    ALIGNED_MALLOC(buf2_mont, mont64, ntt_len);
    ALIGNED_MALLOC(buf3_mont, mont64, ntt_len);
    ALIGNED_MALLOC(tmp_mont, mont64, ntt_len);

    if (buf1_mont == NULL || buf2_mont == NULL || buf3_mont == NULL || tmp_mont == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        abort();
    }

    ntt_short table;
    table.ntt_len = (ntt_len < long_threshold) ? ntt_len : long_threshold;
    table.log_len = log2_64(table.ntt_len);
    //table.omega = g_a;
    //table.iomega = g_b;
    ALIGNED_MALLOC(table.omega, mont64, table.ntt_len);
    ALIGNED_MALLOC(table.iomega, mont64, table.ntt_len);


    for (size_t ii = len1; ii < ntt_len; ii++) {
        buf1_mont[ii] = 0;
        buf2_mont[ii] = 0;
        buf3_mont[ii] = 0;
    }
    for (size_t ii = 0; ii < len1; ii++) {
        buf1_mont[ii] = in1[ii];
        buf2_mont[ii] = in1[ii];
        buf3_mont[ii] = in1[ii];
        _mont_tomont_func(buf1_mont[ii], 1);
        _mont_tomont_func(buf2_mont[ii], 2);
        _mont_tomont_func(buf3_mont[ii], 3);
    }

    for (size_t ii = len2; ii < ntt_len; ii++) {
        tmp_mont[ii] = 0;
    }
    for (size_t ii = 0; ii < len2; ii++) {
        tmp_mont[ii] = in2[ii];
        _mont_tomont_func(tmp_mont[ii], 1);
    }

    clock_t start, end;
    double elapsed;
    cover_nttshort_func(table.log_len, &table, 1);
    
    start = clock();
    convolution_rec_func(buf1_mont, tmp_mont, buf1_mont, &table, ntt_len, 1);
    end = clock();
    elapsed = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
    time_ntt[index][0] = elapsed;

    for (size_t ii = len2; ii < ntt_len; ii++) {
        tmp_mont[ii] = 0;
    }
    for (size_t ii = 0; ii < len2; ii++) {
        tmp_mont[ii] = in2[ii];
        _mont_tomont_func(tmp_mont[ii], 2);
    }
    cover_nttshort_func(table.log_len, &table, 2);

    start = clock();

    convolution_rec_func(buf2_mont, tmp_mont, buf2_mont, &table, ntt_len, 2);
    end = clock();
    elapsed = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
    time_ntt[index][1] = elapsed;

    for (size_t ii = len2; ii < ntt_len; ii++) {
        tmp_mont[ii] = 0;
    }
    for (size_t ii = 0; ii < len2; ii++) {
        tmp_mont[ii] = in2[ii];
        _mont_tomont_func(tmp_mont[ii], 3);
    }
    cover_nttshort_func(table.log_len, &table, 3);

    start = clock();
    convolution_rec_func(buf3_mont, tmp_mont, buf3_mont, &table, ntt_len, 3);
    end = clock();
    elapsed = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
    time_ntt[index][2] = elapsed;

    ALIGNED_FREE(tmp_mont);

    start = clock();
    u192 carry = {0, 0, 0};
    for (size_t ii = 0; ii < conv_len; ii++) {
        u192 temp = {0, 0, 0};
        crt3(buf1_mont[ii], buf2_mont[ii], buf3_mont[ii], temp);
        _u192add(carry, temp);
        out[ii] = (*(carry));
        carry[0] = carry[1];
        carry[1] = carry[2];
        carry[2] = 0;
    }

    out[conv_len] = carry[0];
    end = clock();
    elapsed = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
    time_ntt[index][3] = elapsed;

    ALIGNED_FREE(buf1_mont);
    ALIGNED_FREE(buf2_mont);
    ALIGNED_FREE(buf3_mont);
    ALIGNED_FREE(table.omega);
    ALIGNED_FREE(table.iomega);
}



double test_mul_time(int len1, int len2, int index) {
    u64* in1 = (u64*)malloc(len1 * sizeof(u64));
    u64* in2 = (u64*)malloc(len2 * sizeof(u64));
    u64* out = (u64*)malloc((len1 + len2) * sizeof(u64));
    u64 out_len = len1 + len2;

    srand(1);

    for (int i = 0; i < len1; i++) {
        in1[i] = (u64)rand() * rand();
    }
    for (int i = 0; i < len2; i++) {
        in2[i] = (u64)rand() * rand();
    }

    clock_t start = clock();
    abs_mul64(in1, len1, in2, len2, out, index);
    clock_t end = clock();

    double elapsed = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
    time_ntt[index][4] = elapsed;

    // for (int i = out_len - 10; i < out_len; i++) {
    //     printf("%llu\n", out[i]);
    // }

    free(in1);
    free(in2);
    free(out);
}


void write_csv() {
    FILE* fp = fopen("time_ntt_data.csv", "w");
    if (fp == NULL) {          
        perror("failed to open file"); 
        return;
    }

    // 3. 写入二维数组到CSV文件
    for (int row = 0; row < 251; row++) { 
        for (int col = 0; col < 5; col++) { 
            fprintf(fp, "%.6lf", time_ntt[row][col]);
            if (col != 4) {
                fprintf(fp, ",");
            }
        }
        fprintf(fp, "\n");
    }

    // 4. 关闭文件（必须关闭，确保数据完整写入磁盘）
    fclose(fp);
    printf("successfully write time_ntt_data.csv\n");
}
int main() {

    for (int i = 0; i < 251; i++) {
        int len = 1000 + 39996 * i;
        printf("%d\n", i);
        test_mul_time(len, len, i);
    }
    write_csv();
    return 0;
}
