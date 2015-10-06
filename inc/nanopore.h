#ifndef NANOPORE
#define NANOPORE

typedef struct _nanoporeReadAdjustmentParameters {
    double scale;
    double shift;
    double var;
    double scale_sd;
    double var_sd;
} NanoporeReadAdjustmentParameters;

typedef struct _nanoporeRead {
    int64_t readLength; // 2D read length in nucleotides
    int64_t nbTemplateEvents; // number of events in the array
    int64_t nbComplementEvents; // same
    NanoporeReadAdjustmentParameters templateParams;
    NanoporeReadAdjustmentParameters complementParams;

    char *twoDread; // dread indeed

    int64_t *templateEventMap;
    double *templateEvents; // mean, stdev, length

    int64_t *complementEventMap;
    double *complementEvents;
} NanoporeRead;

NanoporeRead *loadNanoporeReadFromFile(const char *nanoporeReadFile);

#endif