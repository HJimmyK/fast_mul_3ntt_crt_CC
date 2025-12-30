# fast_mul_3ntt_crt_CC

A fast implementation of 3ntt-crt multiplication with C code

The original CPP code can be viewed in this [repository](https://github.com/With-Sky/HintFFT).

This C code is a modified version of the original CPP code, which is optimized for performance and memory usage. The main changes are:

1. Reduced the issue of excessive packaging in Montgomery Domain.
2. Calculate 128 and 192 bit unsigned numbers through macro definition.
3. Explicitly instantiated a large number of template functions through macro definitions

The main performance improvements are:
+ the computational performance of Montgomery domain (approximately 5%)

There is almost no difference in overall performance.

This is because the main bottleneck of the 3ntt-crt algorithm lies in memory read and write operations

The latest measurement results are as follows:
Approximately 10% performance improvement

<img width="1086" height="864" alt="image" src="https://github.com/user-attachments/assets/c6f9dcad-c896-4084-97d0-acc551070b93" />

<img width="1061" height="859" alt="image" src="https://github.com/user-attachments/assets/33346a2f-baba-4048-8bbf-ab4a7c1f83ac" />
</br></br>
The old data: 
<img width="895" height="743" alt="2" src="https://github.com/user-attachments/assets/55895879-170f-4664-b6b2-efae085f856f" />
<img width="1021" height="834" alt="1" src="https://github.com/user-attachments/assets/125cbc41-683c-467f-88c0-67c62b10e86a" />
<img width="983" height="827" alt="3" src="https://github.com/user-attachments/assets/56c9d03b-bb51-4d79-ad18-d38c51a70ba0" />




