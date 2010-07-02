/* This is a truly terrible program and hopefully can go away as soon
  as possible.  Cmake totally sucks at scripting things. When I try to
  do this with regex, cmake consumes the semicolons since it uses that
  for lists. I can't figure out how to work around it without
  documentation / google, and I'm on a train. I would use Haskell (and
  did), (most of what this does is literally 1 line), but a GHC build
  dependence just for that is kind of dump. Really the main problem is
  any of the proper tools for this (or anything else for that matter)
  don't really work on Windows. This is pretty shitty and I wrote it
  partially because I was missing documentation for other things and
  had no internet on the train.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/* I'm a test comment since this program tests itself */

/*/ maybetri//cky shit//*/

typedef enum
{
    REGULAR,
    C_COMMENT,
    CPP_COMMENT
} Type;

typedef struct _Track
{
    size_t start;
    size_t end;
    Type type;
    struct _Track* next;
    struct _Track* prev;
} Track;

void writeInlineFile(const char* bufp, const char* outfile, const char* strName)
{
    FILE* fh;
    FILE* fc;

    size_t len = strlen(outfile) + 1;
    char* hdrName = malloc(len);
    char* cName = malloc(len);

    strcpy(hdrName, outfile);
    strcat(hdrName, ".h");
    strcpy(cName, outfile);
    strcat(cName, ".c");


    fc = fopen(cName, "w");
    if (!fc)
    {
        perror("Open output file");
        exit(EXIT_FAILURE);
    }

    fh = fopen(hdrName, "w");
    if (!fh)
    {
        perror("Open output file");
        exit(EXIT_FAILURE);
    }

    fprintf(fc,
            "\n#include \"%s\"\n"
            "\nconst char* %s = \"%s\""
            "\n\n",
            hdrName,
            strName,
            bufp);

    fprintf(fh,
            "\n#ifndef _%s_h_\n"
            "#define _%s_h_\n\n"
            "extern const char* %s\n\n"
            "#endif /* _%s_h_ */\n\n",
            outfile,
            outfile,
            strName,
            outfile);


    free(hdrName);
    free(cName);
    fclose(fc);
    fclose(fh);
}

const char* showType(Type st)
{
    const char* table[] =
    {
        [REGULAR] = "Regular",
        [C_COMMENT] = "C Commment",
        [CPP_COMMENT] = "C++ Commment"
    };

    if (st > CPP_COMMENT)
        return "Invalid type";

    return table[st];
}

/* Skips over chunks of \0's also */
char* escapeThings(char* buf, size_t size)
{
    char* bufp;
    size_t i, n;

    /* Twice the original size will be big enough, so don't bother
     * figuring out a new one. So lazy. */
    bufp = calloc(sizeof(char), 2 * size);

    for ( i = n = 0; n < size; ++n )
    {
        /* TODO: Other things that need escape */
        if (buf[n] == '\n') /* Escape the newlines */
        {
            bufp[i++] = '\\';
            bufp[i++] = 'n';
        }
        else if (buf[n] == '\r')
        {
            bufp[i++] = '\\';
            bufp[i++] = 'r';
        }
        else if (buf[n] == '\t')
        {
            bufp[i++] = '\\';
            bufp[i++] = 't';
        }
        else if (buf[n] == '\"')
        {
            bufp[i++] = '\\';
            bufp[i++] = '\"';
        }
        else if (buf[n] != '\0')
            bufp[i++] = buf[n];
    }

    return bufp;
}

void eraseComments(Track* cs, char* buf)
{
    Track* c = cs->next; /* Root is junk */

    /* This is a truly terrible program. Erase the comments in the
     * source buffer. We'll then just skip over them after */
    while (c)
    {
        memset(&buf[c->start], 0, c->end - c->start);
        c = c->next;
    }
}

/* Advance to the end of the comment */
char* handleCComment(char* c)
{
    while (*c)
    {
        if (*c == '*')
        {
            ++c;
            if (*c == '/')
                break;
        }
        ++c;
    }
//sample
    return c;
}

/* Advance to the end of the comment */
char* handleCPPComment(char* c)
{
    /* TODO: handle silly windows \r ?? */
    while (*c && *c != '\n')
        ++c;

    return --c;
}

#define SNOC(x,xs)                                               \
        {                                                        \
            (xs)->next = (x);                                    \
            (x)->prev = (xs);                                    \
            (xs) = (x);                                          \
        }                                                        \


void freeTree(Track* ts)
{
    Track* t;
    Track* tmp;

    t = ts;

    while (t)
    {
        tmp = t->next;
        free(t);
        t = tmp;
    }
}


Track* findComments(char* buf)
{
    char* c;
    char* co;

    Track* t;
    Track* ts;
    Track* root = calloc(sizeof(Track), 1);

    c = buf;
    ts = root;

    while (*c)
    {
        switch (*c)
        {
            case '/':
                co = c;
                ++c;
                if (*c == '*')
                {
                    t = calloc(sizeof(Track), 1);
                    t->type = C_COMMENT;
                    t->start = co - buf;
                    c = handleCComment(++c);
                    t->end = c - buf + 1;
                    SNOC(t, ts);
                }
                else if (*c == '/')
                {
                    t = calloc(sizeof(Track), 1);
                    t->type = C_COMMENT;
                    t->start = co - buf;
                    c = handleCPPComment(++c);
                    t->end = c - buf + 1;
                    SNOC(t, ts);
                }
                else
                    --c; /* Backtrack */
                break;
        }
        ++c;
    }

    return root;
}

void arst(const char* filename, const char* outfile, const char* strName)
{
    int size;
    FILE* f;
    char* buf;
    char* bufp;
    Track* cs;

    f = fopen(filename, "r");
    if (f == NULL)
    {
        perror("open input file");
        exit(EXIT_FAILURE);
    }

    fseek(f, 0L, SEEK_END);
    size = ftell(f);
    fseek(f, 0L, SEEK_SET);

    buf = calloc(sizeof(char), size + 2);

    fread(buf, sizeof(char), size + 1, f);

    if (!feof(f))
        perror("Error reading file");

    fclose(f);

    cs = findComments(buf);
    eraseComments(cs, buf);
    freeTree(cs);

    bufp = escapeThings(buf, size);

    writeInlineFile(bufp, outfile, strName);

    free(buf);
    free(bufp);

}


int main(int argc, const char* argv[])
{
    if (argc < 4)
    {
        fprintf(stderr, "Usage: %s: (input file) (base output filename) (string name)\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    arst(argv[1], argv[2], argv[3]);

    return 0;
}
//I'm an edge case
