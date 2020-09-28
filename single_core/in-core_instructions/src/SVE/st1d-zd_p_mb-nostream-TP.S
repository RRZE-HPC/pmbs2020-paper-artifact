#define INSTR st1d
#define NINST 7
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

        mov     x0, sp
        and     x0, x0, #0xFFFFFFFFFFFFFFC0
 
        sub     x1, x0, #512
        sub     x2, x0, #896
        sub     x3, x0, #1024
        sub     x5, x0, #2048
        sub     x6, x0, #4096
        sub     x7, x0, #8192

loop:
        subs      x4, x4, #1
        INSTR    {z0.d}, p0, [x0]
        INSTR    {z1.d}, p0, [x6]
        INSTR    {z2.d}, p0, [x2]
        INSTR    {z3.d}, p0, [x7]
        INSTR    {z4.d}, p0, [x5]  
        INSTR    {z5.d}, p0, [x1]
        INSTR    {z5.d}, p0, [x3]
        bne       loop
done:
        # pop callee-save registers from stack
        ld1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]
        add            sp, sp, #64
        ld1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        add            sp, sp, #64
        
        ret

.size latency, .-latency
