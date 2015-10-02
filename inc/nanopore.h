#ifndef NANOPORE
#define NANOPORE

// any includes?

typedef struct _nanoporeRead {
    int64_t length;
    int64_t nb_events;
    double scale;
    double shift;
    char *read;
    int64_t *eventMap;
    double *events;
} NanoporeRead;

NanoporeRead *loadNanoporeReadFromFile(const char *nanoporeReadFile);

#endif