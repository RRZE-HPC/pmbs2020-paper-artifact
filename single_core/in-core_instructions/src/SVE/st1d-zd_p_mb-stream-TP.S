#define INSTR st1d
#define NINST 8
#define N x0

.arch armv8.2-a+sve
.globl ninst
.data
ninst:
.long NINST
.text
.globl latency
.type latency, @function
.align 2
latency:

        # push callee-save registers onto stack
        sub            sp, sp, #64
        st1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        sub            sp, sp, #64
        st1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]

        mov     x4, N

        fmov    v0.2d, #1.00000000
        fmov    v1.2d, #1.00000000
        fmov    v2.2d, #1.00000000

        ptrue   p0.d
        ptrue   p1.d
        fcpy    z0.d, p0/m, #1.00000000
        fcpy    z1.d, p0/m, #1.00000000
        fcpy    z2.d, p0/m, #1.00000000

        mov     x1, sp
        and     x1, x1, #0xFFFFFFFFFFFFFFC0
        mov     x2, sp
        mov     x3, sp
        mov     x5, sp
        mov     x6, sp
        mov     x7, sp
        sub     x1, x1, #640
        sub     x2, x2, #128
        sub     x3, x3, #192
        sub     x5, x5, #256
        sub     x6, x6, #320
        sub     x7, x7, #2048
loop:
        subs      x4, x4, #1
        st1d	z4.d, p0, [x1]
    	st1d	z5.d, p0, [x1, #1, mul vl]
    	st1d	z6.d, p0, [x1, #2, mul vl]
    	st1d	z7.d, p0, [x1, #3, mul vl]
    	st1d	z12.d, p0, [x1, #4, mul vl]
    	st1d	z13.d, p0, [x1, #5, mul vl]
    	st1d	z14.d, p0, [x1, #6, mul vl]
    	st1d	z15.d, p0, [x1, #7, mul vl]
        bne       loop
done:
        # pop callee-save registers from stack
        ld1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]
        add            sp, sp, #64
        ld1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        add            sp, sp, #64
        
        ret

.size latency, .-latency
