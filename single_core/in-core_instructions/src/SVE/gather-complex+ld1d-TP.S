#define INSTR ld1d
#define NINST 10
#define STRIDE #7
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

        mov     z20.d, p0/z, STRIDE
        mov     z21.d, p0/z, STRIDE
        mov     z22.d, p0/z, STRIDE
        mov     z23.d, p0/z, STRIDE
        mov     z24.d, p0/z, STRIDE
        mov     z25.d, p0/z, STRIDE
        mov     z26.d, p0/z, STRIDE
        mov     z27.d, p0/z, STRIDE
        mov     z28.d, p0/z, STRIDE
        mov     z29.d, p0/z, STRIDE

    
        mov     x0, sp
        # align to 64 byte
        and     x0, x0, #0xFFFFFFFFFFFFFFC0
        sub     x6, x1, #8192
        
        ldr     x1, [sp]
        st1d   {z0.d}, p0, [sp]

loop:
        subs      x4, x4, #1
        INSTR   z0.d, p0/z, [x0, z20.d]
        ld1d    z11.d, p1/z, [x6]
        INSTR   z1.d, p0/z, [x0, z21.d]
        ld1d    z12.d, p1/z, [x6]
        INSTR   z2.d, p0/z, [x0, z22.d]
        ld1d    z13.d, p1/z, [x6]
        INSTR   z3.d, p0/z, [x0, z23.d]
        ld1d    z14.d, p1/z, [x6]
        INSTR   z4.d, p0/z, [x0, z24.d]
        ld1d    z15.d, p1/z, [x6]
        INSTR   z5.d, p0/z, [x0, z25.d] 
        ld1d    z16.d, p1/z, [x6]
        INSTR   z6.d, p0/z, [x0, z26.d]
        ld1d    z17.d, p1/z, [x6]
        INSTR   z7.d, p0/z, [x0, z27.d] 
        ld1d    z18.d, p1/z, [x6]
        INSTR   z8.d, p0/z, [x0, z28.d]
        ld1d    z19.d, p1/z, [x6]
        INSTR   z9.d, p0/z, [x0, z29.d] 
        ld1d    z30.d, p1/z, [x6]
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
