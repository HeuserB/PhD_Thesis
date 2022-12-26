#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>


int main(int argc, char** argv){
	int fd = open(argv[1],O_RDONLY, S_IRUSR | S_IWUSR);
	struct stat sb;

	if (fstat(fd,&sb) == -1)
	{
		perror("couldn't get file size.\n");
	}
	printf("file size is %ld\n", sb.st_size);

	char *file_in_memory =	static_cast<char*>(mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));

	printf("Printing file content, as an array of characters... \n\n");
	for (int i=0; i<sb.st_size; i++)
	{
		printf("%c", file_in_memory[i]);
	}
	printf("\n");
}
