#define INSTR ld1d
#define NINST 6
#define STRIDE #32
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
        fcpy    z6.d, p0/m, #1.00000000

        mov     z0.d, p0/z, STRIDE
        mov     z1.d, p0/z, STRIDE
    
        ldr     x1, [sp]
        st1d   {z0.d}, p0, [sp]

loop:
        subs      x4, x4, #1
        INSTR    z0.d, p0/z, [sp, z1.d]
        INSTR    z1.d, p0/z, [sp, z0.d]
        INSTR    z0.d, p0/z, [sp, z1.d]
        INSTR    z1.d, p0/z, [sp, z0.d]
        INSTR    z0.d, p0/z, [sp, z1.d]
        INSTR    z1.d, p0/z, [sp, z0.d] 
 
        bne       loop
done:
        str     x1, [sp]
        # pop callee-save registers from stack
        ld1            {v12.2d, v13.2d, v14.2d, v15.2d}, [sp]
        add            sp, sp, #64
        ld1            {v8.2d, v9.2d, v10.2d, v11.2d}, [sp]
        add            sp, sp, #64
        
        ret

.size latency, .-latency
