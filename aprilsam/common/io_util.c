/* Copyright (C) 2013-2019, The Regents of The University of Michigan.
All rights reserved.

This software was developed in the APRIL Robotics Lab under the
direction of Edwin Olson, ebolson@umich.edu. This software may be
available under alternative licensing terms; contact the address above.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, see <http://www.gnu.org/licenses/>.
*/

/*$LICENSE*/

/* $COPYRIGHT_UM
$LICENSE_LGPL
*/

#include <assert.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/ioctl.h>
#include <sys/poll.h>
#include <arpa/inet.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

#include "io_util.h"

int ioutils_write_file(const char *path, const void *_buf, long bufsz)
{
    const char *buf = _buf;
    int ret = 0;

    FILE *f = fopen(path, "w");
    if (!f)
        return -1;

    if (bufsz != fwrite(buf, 1, (size_t) bufsz, f))
        ret = -2;

    fclose(f);
    return ret;
}

// returns zero on success. On failure, does not allocate any memory.
int ioutils_read_file(const char *path, void *_buf, long *_bufsz, long maxalloc)
{
    FILE *f = fopen(path, "r");
    if (!f)
        return -1;

    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (maxalloc >= 0 && sz > maxalloc) {
        fclose(f);
        return -2;
    }

    void *buf = malloc(sz + 1);
    if (sz != fread(buf, 1, sz, f))
        assert(0);

    fclose(f);

    ((char*) buf)[sz] = 0;

    *((void**) _buf) = buf;
    *_bufsz = sz;

    return 0;
}

char *
get_unique_filename (const char * path, const char * basename,
                     uint8_t time_prefix, const char * extension)
{
    char prefix[256];

    if (time_prefix) {
        time_t t = time (NULL);
        struct tm ti;
        localtime_r (&t, &ti);

        if (path)
            snprintf (prefix, sizeof (prefix), "%s/%d-%02d-%02d-%s",
                      path, ti.tm_year+1900, ti.tm_mon+1, ti.tm_mday, basename);
        else
            snprintf (prefix, sizeof (prefix), "%d-%02d-%02d-%s",
                      ti.tm_year+1900, ti.tm_mon+1, ti.tm_mday, basename);
    }
    else {
        if (path)
            snprintf (prefix, sizeof (prefix), "%s/%s",
                      path, basename);
        else
            snprintf (prefix, sizeof (prefix), "%s", basename);
    }

    int maxlen = sizeof (prefix) + 4;
    if (extension)
        maxlen += strlen (extension);
    char * filename = malloc (maxlen);

    /* Loop through possible file names until we find one that doesn't already
     * exist.  This way, we never overwrite an existing file. */
    int res;
    int filenum = 0;
    do {
        struct stat statbuf;
        if (extension)
            snprintf (filename, maxlen, "%s.%02d.%s", prefix, filenum,
                      extension);
        else
            snprintf (filename, maxlen, "%s.%02d", prefix, filenum);
        res = stat (filename, &statbuf);
        filenum++;
    } while (res == 0 && filenum < 100);

    if (filenum >= 100) {
        free (filename);
        return NULL;
    }

    if (errno != ENOENT) {
        free (filename);
        return NULL;
    }

    return filename;
}

int fdprintf(int fd, const char *fmt, ...)
{
    int bufsize = 1024;
    char *buf = malloc(bufsize);
    va_list ap;

  again:

    va_start(ap, fmt);
    int n = vsnprintf(buf, bufsize, fmt, ap);
    va_end(ap);

    if (n >= bufsize) {
        bufsize = n + 1;
        buf = realloc(buf, bufsize);
        goto again;
    }

    int res = write_fully(fd, buf, n);

    free(buf);
    return res;
}

int write_fully(int fd, const void *b, int len)
{
    int cnt=0;
    int thiscnt;
    unsigned char *bb=(unsigned char*) b;

    while (cnt<len)
    {
        thiscnt=write(fd, &bb[cnt], len-cnt);
        if (thiscnt<0) {
            return -1;
        }
        cnt+=thiscnt;
    }

    return cnt;
}

int read_fully(int fd, void *b, int len)
{
    int cnt=0;
    int thiscnt;
    unsigned char *bb=(unsigned char*) b;

    while (cnt<len)
    {
        thiscnt=read(fd, &bb[cnt], len-cnt);
        if (thiscnt<0) {
            perror ("read_fully");
            return -1;
        }
        if (thiscnt == 0) {
            fprintf (stderr, "read_fully end of file\n");
            return -1;
        }

        cnt+=thiscnt;
    }

    return cnt;
}

/* returns -1 on error, 0 on timeout, nchar on success. */
int read_timeout(int fd, void *buf, int maxlen, int timeout_ms)
{
    struct pollfd pfd;

    pfd.fd = fd;
    pfd.events = POLLIN;

    int res = poll(&pfd, 1, timeout_ms);
    if (res < 0) // error
        return -1;

    if (res == 0) { // timeout
        errno = ETIMEDOUT;
        return 0;
    }

    ssize_t readlen = read(fd, buf, maxlen);
    if (readlen < 0) {
        return -1;
    }

    if (readlen == 0) {
        errno = EOF;
        return -1;
    }

    return readlen;
}

/* returns -1 on error, nchar on timeout or success. */
int read_fully_timeout(int fd, void *bufin, int len, int timeout_ms)
{
    char *buf = (char*) bufin;
    int readsofar = 0;
    int thisread = 0;

    while (readsofar < len)
    {
        thisread = read_timeout(fd, &buf[readsofar], len-readsofar, timeout_ms);
        if (thisread == 0)
            return 0;
        if (thisread < 0)
            return thisread;

        readsofar += thisread;
    }

    return len;
}

// reads number of bytes available, <0 on error.
int read_available(int fd)
{
    long avail = 0;

    if(ioctl(fd, FIONREAD, &avail) == 0)
        return avail;
    else
        return -1;
}

// read (and discard) all available data.
void read_flush(int fd)
{
    int avail;
    char buf[1024];

    do
    {
        avail=read_available(fd);
        if (avail>1024)
            avail=1024;

        int res = read(fd, buf, avail);
        res++; // trick to ignore compiler warnings about unread result

    } while (avail>0);
}

// returns length of string, -1 on error, 0 on timeout.
int read_line_timeout(int fd, void *buf_in, int maxlen, int timeout_ms)
{
    int timed_out;
    return read_line_timeout_ex(fd, buf_in, maxlen, timeout_ms, &timed_out);
}

int read_line_timeout_ex(int fd, void *buf_in, int maxlen, int timeout_ms, int *timed_out)
{
    int len=0;
    int thislen;
    char *buf=(char*) buf_in;

    *timed_out = 0;

    while (len<maxlen)
    {
        thislen=read_timeout(fd, &buf[len], 1, timeout_ms);
        if (thislen < 0)
            return thislen;
        if (thislen==0) {
            *timed_out = 1;
            return 0;
        }

        if (buf[len]=='\n' || buf[len]=='\r')
        {
            buf[len]=0;
            return len;
        }

        len++;
    }

    buf[len]=0;
    return len;
}

// returns 0 if a line ending in one of the newlinechars was found. *buflen will be set with the total number of characters in buf.
// returns -1 on timeout
// returns -2 on io error
// returns -3 on buffer too short
int read_line_timeout_flags(int fd, char **_buf, int *_buflen, int maxlen, int timeout_ms,
                            const char *newlinechars, const char *stripchars)
{
    int alloc = 16;
    int len = 0;
    char *buf = malloc(alloc);

    while (maxlen < 0 || len < maxlen) {
        char c;
        int thislen = read_timeout(fd, &c, 1, timeout_ms);

        if (thislen == 0) {
            free(buf);
            return -1;
        }

        if (thislen < 0) {
            free(buf);
            return -2;
        }

        int is_stripchar = 0;
        for (const char *p = stripchars; *p; p++)
            is_stripchar |= (c == *p);

        if (!is_stripchar) {
            if (alloc <= len + 1) {
                alloc *= 2;
                buf = realloc(buf, alloc);
            }
            buf[len++] = c;
        }

        int is_newline = 0;
        for (const char *p = newlinechars; *p; p++)
            is_newline |= (c == *p);

        if (is_newline) {
            buf[len] = 0;
            *_buf = buf;
            *_buflen = len;
            return 0;
        }
    }

    free(buf);
    return -3;
}

int fwrite64(FILE *f, int64_t v64)
{
    int32_t v = v64>>32;
    if (0 != fwrite32(f, v)) return -1;
    v = v64 & 0xffffffff;
    return fwrite32(f, v);
}

int fwrite32(FILE *f, int32_t v)
{
    v = htonl(v);
    if (fwrite(&v, 4, 1, f) == 1) return 0; else return -1;
}

int fread32(FILE *f, int32_t *v32)
{
    int32_t v;

    if (fread(&v, 4, 1, f) != 1)
        return -1;

    *v32 = ntohl(v);

    return 0;
}

int fread64(FILE *f, int64_t *v64)
{
    int32_t v1, v2;

    if (fread32(f, &v1))
        return -1;

    if (fread32(f, &v2))
        return -1;

    *v64 =     (((int64_t) v1)<<32) | (((int64_t) v2)&0xffffffff);

    return 0;
}
