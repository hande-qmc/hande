#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

void alloc_shared_posix(void **ptr, const char *handle, const size_t len)
{
	int fd;

#if defined PARALLEL && defined ENABLE_SHMEM_POSIX
	fd = shm_open(handle, O_CREAT | O_RDWR, S_IRWXU);

	if(fd == -1 ) { perror( "shm_open" ); exit(0);}

    if (ftruncate(fd, len) == -1) { perror("Setting size of shared memory region failed."); }

    // This block of memory *must* be unmapped before the program terminates,
    // otherwise it will remain unusable until the computer is rebooted.
	*ptr = (int*) mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0 );

    if (*ptr == MAP_FAILED) { perror("Mapping of shared memory failed."); }

    /* the mapped region will remain accessible even once the file handle is closed. */
    close(fd);
#endif

}

void free_shared_posix(void *ptr, char *handle, int len)
{
#if defined PARALLEL && defined ENABLE_SHMEM_POSIX
	munmap(ptr, len);
	shm_unlink(handle);
#endif
}
