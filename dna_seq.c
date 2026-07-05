#include <stdio.h>
#include <string.h>

#define MAX_LINE 1000

int main(int argc, char *argv[]) {

    if (argc < 2) {
        printf("Usage: ./dna_analyzer sequence.fasta\n");
        return 1;
    }

    FILE *fp = fopen(argv[1], "r");

    if (fp == NULL) {
        printf("Error: Could not open file '%s'\n", argv[1]);
        return 1;
    }

    char line[MAX_LINE];
    int count_A = 0;
    int count_T = 0;
    int count_G = 0;
    int count_C = 0;

    printf("Reading FASTA file: %s\n\n", argv[1]);

    while (fgets(line, MAX_LINE, fp) != NULL) {

        if (line[0] == '>') {
            printf("Header found: %s", line);
        } else {
    printf("Sequence: %s", line);
            for (int i = 0; line[i] != '\0'; i++) {
                if (line[i] == 'A') count_A++;
                if (line[i] == 'T') count_T++;
                if (line[i] == 'G') count_G++;
                if (line[i] == 'C') count_C++;
            }
        }
    }

    fclose(fp);

    printf("\n--- Nucleotide Count ---\n");
    printf("A: %d\n", count_A);
    printf("T: %d\n", count_T);
    printf("G: %d\n", count_G);
    printf("C: %d\n", count_C);
    printf("Total: %d\n", count_A + count_T + count_G + count_C);

    return 0;
}