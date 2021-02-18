/*
 * =======================================================================================
 *
 *      Author:   Jan Eitzinger (je), jan.eitzinger@fau.de
 *      Copyright (c) 2019 RRZE, University Erlangen-Nuremberg
 *
 *      Permission is hereby granted, free of charge, to any person obtaining a copy
 *      of this software and associated documentation files (the "Software"), to deal
 *      in the Software without restriction, including without limitation the rights
 *      to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *      copies of the Software, and to permit persons to whom the Software is
 *      furnished to do so, subject to the following conditions:
 *
 *      The above copyright notice and this permission notice shall be included in all
 *      copies or substantial portions of the Software.
 *
 *      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *      IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *      FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *      AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *      LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *      OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *      SOFTWARE.
 *
 * =======================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/mman.h>

void* hp_allocate(size_t nbytes) {
    int page_size = 2097152;
    void* ret_ptr = NULL;
    size_t num_large_pages = nbytes / page_size;
    if (nbytes > num_large_pages * page_size) {
        num_large_pages++;
    }
    nbytes = (size_t) num_large_pages * page_size;
    //printf("trying to allocate %ld pages\n", num_large_pages);
    ret_ptr = mmap(NULL, nbytes,
                   PROT_READ | PROT_WRITE,
                   MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB,
                   -1, 0); 
    if ((ret_ptr == (void *)(-1))) {
        fprintf(stderr,"mmap call failed\n");
        exit(1);
    }
    return ret_ptr;
}

void hp_free(void * ptr, size_t nbytes) {
    munmap(ptr, nbytes);
}

void* allocate (int alignment, size_t bytesize)
{
    int errorCode;
    void* ptr;

    errorCode =  posix_memalign(&ptr, alignment, bytesize);

    if (errorCode) {
        if (errorCode == EINVAL) {
            fprintf(stderr,
                    "Error: Alignment parameter is not a power of two\n");
            exit(EXIT_FAILURE);
        }
        if (errorCode == ENOMEM) {
            fprintf(stderr,
                    "Error: Insufficient memory to fulfill the request\n");
            exit(EXIT_FAILURE);
        }
    }

    if (ptr == NULL) {
        fprintf(stderr, "Error: posix_memalign failed!\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}
