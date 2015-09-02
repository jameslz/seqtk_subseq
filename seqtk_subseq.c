#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khash.h"
KHASH_SET_INIT_STR(reg)

#include "kstring.h"

void print_help(const char *command_line) {
    printf("Usage: %s <fasta|fastq> <list> <flag: 0|1 >\n", command_line);
}

inline void print_seq(const kseq_t *s){

    if( s->qual.l ){
        printf("@%s\n%s\n+\n%s\n", s->name.s, s->seq.s, s->qual.s);
    }else{
        printf(">%s", s->name.s);
        if (s->comment.l) printf(" %s", s->comment.s);
        printf("\n%s\n", s->seq.s);   
    }
}

int main(int argc, char const *argv[]){
    
    if(argc == 1) {
        print_help( argv[0] );
        exit(1);
    }

    khash_t(reg) *h;
    khint_t k;

    int i, absent;

    h = kh_init(reg);

    gzFile     fp; 
    kstream_t  *ks;    
    kstring_t  str  = {0, 0, 0};
    fp = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r"); 
    ks              = ks_init(fp);

    while( ks_getuntil( ks, '\n', &str, 0) >=  0){
        if( str.l == 0 ) continue;
        k = kh_put(reg, h, str.s, &absent);
        if(absent) kh_key(h, k) = strdup(str.s);
    }

    gzclose(fp);

    kseq_t   *seq;
    fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r"); 
    seq      = kseq_init(fp);
    int      l, flag = atoi(argv[3]);

    while ((l = kseq_read(seq)) >= 0) {
    
        k = kh_get(reg, h, seq->name.s);
        
        if(k != kh_end(h) && flag == 1){
            print_seq(seq);
        }else if(k == kh_end(h) && flag == 0){
            print_seq(seq);
        }
    }
    
    kseq_destroy(seq);
    gzclose(fp);

    for (k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free((char*)kh_key(h, k));
        }
    }
    kh_destroy(reg, h);
    return 0;
}
