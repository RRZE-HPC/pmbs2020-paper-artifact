#define INSTR ld1d
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
        fcpy    z0.d, p0/m, #1.00000000
        fcpy    z1.d, p0/m, #1.00000000
        fcpy    z2.d, p0/m, #1.00000000

        mov     x1, sp
        # align to 64 byte
        and     x1, x1, #0xFFFFFFFFFFFFFFC0

        sub     x1, x1, #64 
        sub     x2, x1, #128
        sub     x3, x1, #192
        sub     x5, x1, #256
        sub     x6, x1, #320
        sub     x7, x1, #384
        sub     x8, x1, #448
        sub     x9, x1, #512
loop:
        subs      x4, x4, #1
        INSTR    {z0.d}, p0/z, [x1]
        INSTR    {z1.d}, p0/z, [x2]
        INSTR    {z2.d}, p0/z, [x3]
        INSTR    {z3.d}, p0/z, [x5]
        INSTR    {z4.d}, p0/z, [x6]
        INSTR    {z5.d}, p0/z, [x7] 
        INSTR    {z6.d}, p0/z, [x8] 
        INSTR    {z7.d}, p0/z, [x9] 
        bne       loop
done:
        # pop callee-save registers from stack
        ld1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]
        add            sp, sp, #64
        ld1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        add            sp, sp, #64
        
        ret

.size latency, .-latency
