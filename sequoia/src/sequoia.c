/*----------------------------------------------------------------------
  File    : sequoia.c
  Contents: find frequent sequences with unique item occurrences
  Author  : Christian Borgelt
  History : 2008.08.11 file created from swog.c
            2010.08.12 redesigned to work mainly with arrays
            2010.08.13 alternative variants of closed() added
            2010.08.14 variant without item weights added
            2010.08.22 adapted to modified modules tabread and tract
            2010.11.24 adapted to modified error reporting (tract)
            2010.12.11 adapted to a generic error reporting function
            2011.03.20 optional integer transaction weights added
            2011.07.08 adapted to modified function tbg_recode()
            2011.08.28 output of sequence counters per size added
            2012.07.23 bug in item base creation fixed (IB_WEIGHTS)
            2013.04.04 adapted to type changes in module tract
            2013.10.15 checks of result of isr_isetx() added
            2013.10.18 optional pattern spectrum collection added
            2014.05.12 option -F# added (support border for filtering)
            2014.08.26 adapted to modified item set reporter interface
            2014.10.24 changed from LGPL license to MIT license
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifndef TA_READ
#define TA_READ
#endif
#include "tract.h"
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#include "report.h"
#include "error.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "sequoia"
#define DESCRIPTION "sequence mining with unique occurrences " \
                    "of items and weight averaging"
#define VERSION     "version 2.16 (2016.10.15)        " \
                    "(c) 2010-2016   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid sequence length */
#define E_SUPPORT   (-11)       /* invalid minimum sequence support */
#define E_MEASURE   (-13)       /* invalid evaluation measure */
/* error codes -16 to -25 defined in tract.h */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define CLOCK(t)    ((t) = clock())
#else                           /* if quiet version, */
#define MSG(...)    ((void)0)   /* suppress messages */
#define CLOCK(t)    ((void)0)
#endif

#define SEC_SINCE(t)  ((double)(clock()-(t)) /(double)CLOCKS_PER_SEC)

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- pattern occurrence --- */
  SUPP       wgt;               /* weight of containing transaction */
  const ITEM *items;            /* items  of containing transaction */
  const ITEM **ips;             /* (positions of) items in pattern */
} PATOCC;                       /* (pattern occurrence) */

typedef struct {                /* --- occurrence extension --- */
  const ITEM *item;             /* extension item in transaction */
  PATOCC     *occ;              /* pattern occurrence to extend */
} OCCEXT;                       /* occurrence extension */

typedef struct {                /* --- pattern extension --- */
  SUPP       supp;              /* support of extension item */
  TID        cnt;               /* number of occurrences */
  OCCEXT     *oxs;              /* occurrence extensions */
} PATEXT;                       /* (pattern extension) */

typedef struct {                /* --- pattern occurrence --- */
  SUPP       wgt;               /* weight of containing transaction */
  WITEM      *items;            /* items  of containing transaction */
  WITEM      **ips;             /* (positions of) items in pattern */
} WPATOCC;                      /* (pattern occurrence) */

typedef struct {                /* --- occurrence extension --- */
  WITEM      *item;             /* extension item in transaction */
  WPATOCC    *occ;              /* pattern occurrence to extend */
} WOCCEXT;                      /* occurrence extension */

typedef struct {                /* --- pattern extension --- */
  SUPP       supp;              /* support of extension item */
  TID        cnt;               /* number of occurrences */
  WOCCEXT    *oxs;              /* occurrence extensions */
} WPATEXT;                      /* (pattern extension) */

typedef struct {                /* --- recursion data --- */
  int        target;            /* target type (e.g. closed/maximal) */
  int        mode;              /* operation mode (e.g. pruning) */
  ITEM       cnt;               /* total number of items */
  ITEM       zmax;              /* maximum length of a sequence */
  SUPP       smin;              /* minimum support of a sequence */
  TID        *frqs;             /* item counters (closedness check) */
  ITEM       *buf;              /* item buffer   (closedness check) */
  ITEM       *items;            /* current pattern sequence: items */
  double     *wgts;             /* current pattern sequence: weights */
  ISREPORT   *report;           /* item set/sequence reporter */
} RECDATA;                      /* (recursion data) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined SOIA_MAIN
/* --- error messages --- */
static const char *errmsgs[] = {
  /* E_NONE      0 */  "no error",
  /* E_NOMEM    -1 */  "not enough memory",
  /* E_FOPEN    -2 */  "cannot open file %s",
  /* E_FREAD    -3 */  "read error on file %s",
  /* E_FWRITE   -4 */  "write error on file %s",
  /* E_STDIN    -5 */  "double assignment of standard input",
  /* E_OPTION   -6 */  "unknown option -%c",
  /* E_OPTARG   -7 */  "missing option argument",
  /* E_ARGCNT   -8 */  "wrong number of arguments",
  /* E_TARGET   -9 */  "invalid target type '%c'",
  /* E_SIZE    -10 */  "invalid sequence length %d",
  /* E_SUPPORT -11 */  "invalid minimum support %g",
  /*           -12 */  NULL,
  /* E_MEASURE -13 */  "invalid evaluation measure '%c'",
  /*           -14 */  NULL,
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /*           -16 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef SOIA_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set/sequence reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
static double   *border = NULL; /* support border for filtering */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions for Debugging
----------------------------------------------------------------------*/
#ifndef NDEBUG

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show (PATEXT *exts, ITEM n, ITEM len, int ind)
{                               /* --- show pattern extensions */
  ITEM   i, m;                  /* loop variables */
  TID    k;                     /* loop variable */
  PATOCC *o;                    /* to traverse the pattern occs. */
  ITEM   *s;                    /* to traverse the items */

  assert(exts);                 /* check the function arguments */
  for (i = 0; i < n; i++) {     /* traverse the pattern extensions */
    if (exts[i].cnt <= 0) continue;
    indent(ind);                /* print extension item information */
    printf("%02"ITEM_FMT":%s: ", i, ib_name(ibase, i));
    printf("%"TID_FMT"/%"SUPP_FMT"\n", exts[i].cnt, exts[i].supp);
    for (k = 0; k < exts[i].cnt; k++) {
      o = exts[i].oxs[k].occ;   /* traverse the pattern occurrences */
      indent(ind); printf("  ");/* indent the output line */
      for (m = 0; m < len; m++) {
        s = o->ips[m];          /* traverse the pattern */
        printf(" %s", ib_name(ibase, *s));
      }                         /* print the pattern items */
      printf(" |");             /* print a tail separator */
      for (s = exts[i].oxs[k].item; *s >= 0; s++)
        printf(" %s", ib_name(ibase, *s));
      printf("\n");             /* print the tail items */
    }                           /* and terminate the output line */
  }
}  /* show() */

/*--------------------------------------------------------------------*/

static void xshow (WPATEXT *exts, ITEM n, ITEM len, int ind)
{                               /* --- show pattern extensions */
  int     i, k, m;              /* loop variables */
  WPATOCC *o;                   /* to traverse the pattern occs. */
  WITEM   *x;                   /* to traverse the (extended) items */

  assert(exts);                 /* check the function arguments */
  for (i = 0; i < n; i++) {     /* traverse the pattern extensions */
    if (exts[i].cnt <= 0) continue;
    indent(ind);                /* print extension item information */
    printf("%02"ITEM_FMT":%s: ", i, ib_name(ibase, i));
    printf("%"TID_FMT"/%"SUPP_FMT"\n", exts[i].cnt, exts[i].supp);
    for (k = 0; k < exts[i].cnt; k++) {
      o = exts[i].oxs[k].occ;   /* traverse the pattern occurrences */
      indent(ind); printf("  ");/* indent the output line */
      for (m = 0; m < len; m++) {
        x = o->ips[m];          /* traverse the pattern */
        printf(" %s:%g", ib_name(ibase, x->item), x->wgt);
      }                         /* print the pattern items */
      printf(" |");             /* print a tail separator */
      for (x = exts[i].oxs[k].item; x->item >= 0; x++)
        printf(" %s:%g", ib_name(ibase, x->item), x->wgt);
      printf("\n");             /* print the tail items */
    }                           /* and terminate the output line */
  }
}  /* xshow() */

#endif
/*----------------------------------------------------------------------
  Sequence Mining with Unique Item Occurrences (no item weights)
----------------------------------------------------------------------*/

static int closed (PATEXT *ext, ITEM n, RECDATA *rd)
{                               /* --- check for a closed extension */
  TID        i, k, c;           /* loop variables, buffer */
  const ITEM *s, *z;            /* to traverse the items */
  ITEM       *b;                /* to traverse copied items */
  PATOCC     *o;                /* to traverse pattern occurrences */
  OCCEXT     *x;                /* to traverse occurrence extensions */

  assert(ext                    /* check the function arguments */
  &&    (ext->cnt > 0) && (n > 0) && rd);
  for (k = 0; k < ext->cnt; k++) {   /* add item to the occurrences */
    x = ext->oxs +k; x->occ->ips[n-1] = x->item; }
  b = rd->buf; k = 0;           /* get a buffer for occurring items */
  while (--n >= 0) {            /* traverse the current pattern */
    for (i = 0; i < ext->cnt; i++) {
      o = ext->oxs[i].occ;      /* traverse the pattern occurrences */
      s = (n > 0) ? o->ips[n-1]+1 : o->items;
      z = o->ips[n];            /* get the bounds of the current gap */
      for (k = 0; s < z; s++) { /* traverse the current gap */
        c = ++rd->frqs[*s];     /* count gap items that occur */
        if (c >  i) k++;        /* in all pattern occurrences */
        if (c <= 1) *b++ = *s;  /* collect all occurring items */
      }                         /* (for later counter clearing) */
      if (k <= 0) break;        /* if the item counter is zero, */
    }                           /* no item is in all occurrences */
    while (b > rd->buf)         /* clear the used occurrence counters */
      rd->frqs[*--b] = 0;       /* (faster than always clearing all) */
    if (k > 0) return 0;        /* if an item is in all occurrences, */
  }                             /* the extension is not closed */
  return -1;                    /* otherwise the extension is closed */
}  /* closed() */

/*--------------------------------------------------------------------*/

static SUPP recurse (PATEXT *exts, size_t z, ITEM len, RECDATA *rd)
{                               /* --- recursive pattern search */
  ITEM       i, k;              /* loop variables */
  SUPP       s, max;            /* (maximum) extension support */
  PATEXT     *cond = NULL;      /* conditional pattern extensions */
  PATEXT     *e, *c;            /* to traverse pattern extensions */
  PATOCC     *o;                /* to traverse pattern occurrences */
  OCCEXT     *x;                /* to traverse occurrence extensions */
  const ITEM *p;                /* to traverse the tail items */

  assert(exts                   /* check the function arguments */
  &&    (z > 0) && (len >= 0) && rd);
  if (++len <= rd->zmax) {      /* if the pattern can be extended */
    cond = (PATEXT*)malloc((size_t)rd->cnt *sizeof(PATEXT)
                          +        z       *sizeof(OCCEXT));
    if (!cond) return -1;       /* allocate memory for the pattern */
    x = (OCCEXT*)(cond+rd->cnt);/* and occurrence extensions */
    for (k = 0; k < rd->cnt; k++) { cond[k].oxs = x; x += exts[k].cnt; }
  }                             /* organize the occ. extension arrays */
  for (max = s = 0, i = 0; i < rd->cnt; i++) {
    e = exts +i;                /* traverse the pattern extensions */
    if (e->supp < rd->smin)     /* if extension item is infrequent, */
      continue;                 /* the item need not be processed */
    if (e->supp > max)          /* find maximal extension support */
      max = e->supp;            /* (for test if a pattern is closed) */
    if ((rd->mode & ISR_CLOSED) /* if to find only closed sequences */
    &&  !closed(e, len, rd))    /* and the extension is not closed, */
      continue;                 /* the item need not be processed */
    isr_add(rd->report, i, e->supp);
    if (!cond) s = 0;           /* add current item to the reporter */
    else {                      /* if to compute cond. extensions */
      for (k = 0; k < rd->cnt; k++)   /* clear the item variables */
        cond[k].supp = 0, cond[k].cnt = 0;
      for (z = 0, k = 0; k < e->cnt; k++) {
        x = e->oxs +k;          /* traverse the occurrence extensions */
        o = x->occ;             /* get corresp. pattern occurrence */
        for (p = x->item; *++p >= 0; z++) {
          c = cond   +*p;       /* traverse the tail of the sequence */
          x = c->oxs +c->cnt++; /* append an occurrence extension */
          x->item  = p;         /* to the array for the tail item, */
          x->occ   = o;         /* note the extension item position, */
          c->supp += o->wgt;    /* copy the pattern occurrence and */
        }                       /* sum sequences weights (support) */
      }
      if (z > 0) {              /* if cond. extensions are not empty */
        s = recurse(cond, z, len, rd);
        if (s < 0) break;       /* find frequent patterns recursively */
      }                         /* and check for a recursion error */
    }
    if ((!(rd->mode & ISR_CLOSED)  /* if to report all patterns */
    ||  (s < e->supp))             /* or the pattern is closed */
    && (isr_report(rd->report) < 0)) {
      s = -1; break; }          /* report the current pattern */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  if (cond) free(cond);         /* delete the conditional extensions */
  return (s < 0) ? s : max;     /* return maximal extension support */
}  /* recurse() */

/*--------------------------------------------------------------------*/

static int sequoia (TABAG *tabag, int target, SUPP smin, int mode,
                    ISREPORT *report)
{                               /* --- search for frequent sequences */
  ITEM       i, k;              /* loop variable, number of items */
  TID        j, n;              /* loop variable, number of trans. */
  size_t     z;                 /* number of item instances */
  SUPP       r;                 /* result of recursion */
  TRACT      *t;                /* to traverse the transactions */
  const ITEM *s, **p;           /* to traverse the items */
  OCCEXT     *x;                /* to traverse occurrence extensions */
  PATOCC     *occs, *o;         /* array of pattern occurrences */
  PATEXT     *exts, *e;         /* array of pattern extensions */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.target = target;           /* store target type and search mode */
  rd.mode   = mode;             /* check and adapt minimum support */
  rd.smin   = (smin > 0) ? smin : 1;
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  rd.report = report;           /* initialize the recursion data */
  rd.zmax   = isr_zmax(report); /* (reporter and max. seq. length) */
  rd.cnt    = k = tbg_itemcnt(tabag);   /* get the number of items */
  if (k <= 0) return isr_report(report);
  rd.buf  = (ITEM*)malloc((size_t)k *sizeof(ITEM)
                         +(size_t)k *sizeof(TID));
  if (!rd.buf) return -1;       /* create an item buffer */
  rd.frqs = (TID*)(rd.buf +k);  /* and an item counter array */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  z = tbg_extent(tabag);        /* and the number of item instances */
  occs = (PATOCC*)malloc((size_t)n *sizeof(PATOCC)
                        +(size_t)z *sizeof(ITEM*));
  if (!occs) { free(rd.buf); return -1; }
  p = (const ITEM**)(occs +n);  /* create pattern occurrences */
  memset(rd.frqs, 0, (size_t)k *sizeof(TID));
  for (j = 0; j < n; j++) {     /* traverse the transactions and */
    t = tbg_tract(tabag, j);    /* create a pattern occurrence */
    o = occs +j;                /* for each transaction */
    o->wgt = ta_wgt(t);         /* note the transaction weight and */
    o->ips = p; p += ta_size(t);/* organize extension item arrays */
    for (s = o->items = ta_items(t); *s >= 0; s++)
      rd.frqs[*s]++;            /* note the item array and */
  }                             /* count the item occurrences */
  exts = (PATEXT*)malloc((size_t)k *sizeof(PATEXT)
                        +(size_t)z *sizeof(OCCEXT));
  if (!exts) return -1;         /* allocate memory for pattern */
  x = (OCCEXT*)(exts +k);       /* and occurrence extensions */
  for (i = 0; i < k; i++) {     /* initialize the pattern extensions */
    e = exts+i; e->supp = 0; e->cnt = 0; e->oxs = x; x += rd.frqs[i]; }
  for (j = 0; j < n; j++) {     /* traverse the transactions and */
    o = occs +j;                /* the items in each transaction */
    for (s = o->items; *s >= 0; s++) {
      e = exts +*s;             /* get the corresp. pattern extension */
      x = e->oxs +e->cnt++;     /* and its next occurrence extension */
      x->item  = s;             /* set the extension item */
      x->occ   = o;             /* and the pattern occurrence */
      e->supp += o->wgt;        /* sum transaction weights (support) */
    }                           /* (exts represents the possible */
  }                             /* extensions of the empty sequence) */
  memset(rd.frqs, 0, (size_t)k *sizeof(TID));
  r = recurse(exts, z, 0, &rd); /* search for frequent sequences */
  if ( (r >= 0)                 /* if no error occurred */
  &&  ((r < tbg_wgt(tabag))     /* if the empty sequence is closed */
  ||  !(mode & ISR_CLOSED)))    /* or all sequences are requested, */
    r = isr_report(report);     /* report the empty sequence */
  free(rd.buf);                 /* deallocate the pattern buffer, */
  free(exts); free(occs);       /* pattern extensions and occurrences */
  return (r < 0) ? (int)r : 0;  /* return the error status */
}  /* sequoia() */

/*----------------------------------------------------------------------
  Sequence Mining with Unique Item Occurrences and Weight Averaging
----------------------------------------------------------------------*/

static int closed_iw (WPATEXT *ext, ITEM n, RECDATA *rd)
{                               /* --- check for a closed extension */
  TID     i, k, c;              /* loop variables, buffer */
  ITEM    *b;                   /* to traverse the item buffer */
  WPATOCC *o;                   /* to traverse pattern occurrences */
  WITEM   *x, *z;               /* to traverse the (extended) items */

  assert(ext                    /* check the function arguments */
  &&    (ext->cnt > 0) && (len > 0) && rd);
  b = rd->buf; k = 0;           /* get a buffer for occurring items */
  while (--n >= 0) {            /* traverse the current pattern */
    for (i = 0; i < ext->cnt; i++) {
      o = ext->oxs[i].occ;      /* traverse the pattern occurrences */
      x = (n > 0) ? o->ips[n-1]+1 : o->items;
      z = o->ips[n];            /* get the bounds of the current gap */
      for (k = 0; x < z; x++) { /* traverse the current gap */
        c = ++rd->frqs[x->item];/* count gap items that occur */
        if (c >  i) k++;        /* in all pattern occurrences */
        if (c <= 1) *b++ = x->item;
      }                         /* collect all occurring items */
      if (k <= 0) break;        /* if the item counter is zero, */
    }                           /* no item is in all occurrences */
    while (b > rd->buf)         /* clear the used occurrence counters */
      rd->frqs[*--b] = 0;       /* (faster than always clearing all) */
    if (k > 0) return 0;        /* if an item is in all occurrences, */
  }                             /* the extension is not closed */
  return -1;                    /* otherwise the extension is closed */
}  /* closed_iw() */

/*--------------------------------------------------------------------*/

static SUPP rec_iw (WPATEXT *exts, size_t z, ITEM len, RECDATA *rd)
{                               /* --- recursive pattern search */
  ITEM    i, k, m;              /* loop variables */
  SUPP    s, max;               /* (maximum) extension support */
  WPATEXT *cond = NULL;         /* conditional pattern extensions */
  WPATEXT *e, *c;               /* to traverse pattern extensions */
  WPATOCC *o;                   /* to traverse pattern occurrences */
  WOCCEXT *x;                   /* to traverse occurrence extensions */
  WITEM   *p;                   /* to traverse the tail items */

  assert(exts                   /* check the function arguments */
  &&    (z > 0) && (len >= 0) && rd);
  if (++len <= rd->zmax) {      /* if the pattern can be extended */
    cond = (WPATEXT*)malloc((size_t)rd->cnt *sizeof(WPATEXT)
                           +        z       *sizeof(WOCCEXT));
    if (!cond) return -1;       /* allocate memory for the pattern */
    x = (WOCCEXT*)(cond +rd->cnt);    /* and occurrence extensions */
    for (k = 0; k < rd->cnt; k++) { cond[k].oxs = x; x += exts[k].cnt; }
  }                             /* organize the occ. extension arrays */
  for (max = s = 0, i = 0; i < rd->cnt; i++) {
    e = exts +i;                /* traverse the pattern extensions */
    if (e->supp < rd->smin)     /* if extension item is infrequent, */
      continue;                 /* the item need not be processed */
    if (e->supp > max)          /* find maximal extension support */
      max = e->supp;            /* (for test if a pattern is closed) */
    rd->items[len-1] = i;       /* add the ext. item to the pattern */
    for (k = 0; k < e->cnt; k++) {        /* and to its occurrences */
      x = e->oxs+k; x->occ->ips[len-1] = x->item; }
    if ((rd->mode & ISR_CLOSED) /* if to find only closed sequences */
    &&  !closed_iw(e, len, rd)) /* and the extension is not closed, */
      continue;                 /* the item need not be processed */
    if (!cond) s = 0;           /* if no extensions, clear support */
    else {                      /* if to compute cond. extensions */
      for (k = 0; k < rd->cnt; k++)   /* clear the item variables */
        cond[k].supp = 0, cond[k].cnt = 0;
      for (z = 0, k = 0; k < e->cnt; k++) {
        x = e->oxs +k;          /* traverse the occurrence extensions */
        o = x->occ;             /* get corresp. pattern occurrence */
        for (p = x->item; (++p)->item >= 0; z++) {
          c = cond   +p->item;  /* traverse the tail of the sequence */
          x = c->oxs +c->cnt++; /* append an occurrence extension */
          x->item  = p;         /* to the array for the tail item, */
          x->occ   = o;         /* note the extension item position, */
          c->supp += o->wgt;    /* copy the pattern occurrence and */
        }                       /* sum sequences weights (support) */
      }
      if (z > 0) {              /* if cond. extensions are not empty */
        s = rec_iw(cond, z, len, rd);
        if (s < 0) break;       /* find frequent patterns recursively */
      }                         /* and check for a recursion error */
    }
    if ((rd->mode & ISR_CLOSED) /* if to report only closed patterns */
    &&  (s >= e->supp))         /* and the pattern is not closed, */
      continue;                 /* continue with the next item */
    for (k = 0; k < len; k++)   /* traverse the current pattern and */
      rd->wgts[k] = 0;          /* clear (conditional) item weights */
    for (k = 0; k < e->cnt; k++) {
      o = e->oxs[k].occ;        /* traverse the pattern occurrences */
      for (m = 0; m < len; m++) /* and their item occurrences and */
        rd->wgts[m] += (double)o->wgt *o->ips[m]->wgt;
    }                           /* sum (conditional) item weights */
    if (isr_isetx(rd->report,rd->items,len,rd->wgts,e->supp,0,0) < 0) {
      s = -1; break; }          /* report the current pattern */
  }
  if (cond) free(cond);         /* delete the conditional extensions */
  return (s < 0) ? s : max;     /* return the error status or */
}  /* rec_iw() */               /* the maximal extension support */

/*--------------------------------------------------------------------*/

static int sequoia_iw (TABAG *tabag, int target, SUPP smin,
                       int mode, ISREPORT *report)
{                               /* --- search for frequent sequences */
  ITEM    i, k;                 /* loop variable, number of items */
  TID     j, n;                 /* loop variable, number of trans. */
  size_t  z;                    /* number of item instances */
  SUPP    r;                    /* result of recursion */
  WTRACT  *t;                   /* to traverse the transactions */
  WITEM   *s, **p;              /* to traverse the (extended) items */
  WOCCEXT *x;                   /* to traverse occurrence extensions */
  WPATOCC *occs, *o;            /* array of pattern occurrences */
  WPATEXT *exts, *e;            /* array of pattern extensions */
  RECDATA rd;                   /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.target = target;           /* store target type and search mode */
  rd.mode   = mode;             /* check and adapt minimum support */
  rd.smin   = (smin > 0) ? smin : 1;
  if (tbg_wgt(tabag) < rd.smin) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  rd.report = report;           /* initialize the recursion data */
  rd.zmax = isr_zmax(report);   /* (reporter and max. seq. length) */
  rd.cnt  = k = tbg_itemcnt(tabag);
  if (k <= 0)                   /* get and check the number of items */
    return (isr_isetx(report, NULL, 0, NULL, tbg_wgt(tabag), 0, 0) < 0)
         ? -1 : 0;              /* report the empty sequence */
  rd.wgts = (double*)malloc((size_t) k    *sizeof(double)
                           +(size_t)(k+k) *sizeof(ITEM)
                           +(size_t) k    *sizeof(TID));
  if (!rd.wgts) return -1;      /* create a pattern weight array and */
  rd.items = (ITEM*)(rd.wgts+k);/* an item buffer for the reporting, */
  rd.buf   = rd.items +k;       /* a buffer for the closedness check */
  rd.frqs  = (TID*) (rd.buf +k);/* and an item counter array */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  z = tbg_extent(tabag);        /* and the number of item instances */
  occs = (WPATOCC*)malloc((size_t)n *sizeof(WPATOCC)
                         +(size_t)z *sizeof(WITEM*));
  if (!occs) { free(rd.wgts); return -1; }
  p = (WITEM**)(occs +n);       /* create pattern occurrences */
  memset(rd.frqs, 0, (size_t)k *sizeof(TID));
  for (j = 0; j < n; j++) {     /* traverse the transactions and */
    t = tbg_wtract(tabag, j);   /* create a pattern occurrence */
    o = occs +j;                /* for each transaction */
    o->wgt = wta_wgt(t);        /* note the transaction weight and */
    o->ips = p; p += wta_size(t);/* organize extension item arrays */
    for (s = o->items = wta_items(t); s->item >= 0; s++)
      rd.frqs[s->item]++;       /* note the item array and */
  }                             /* count the item occurrences */
  exts = (WPATEXT*)malloc((size_t)k *sizeof(WPATEXT)
                         +(size_t)z *sizeof(WOCCEXT));
  if (!exts) return -1;         /* allocate memory for pattern */
  x = (WOCCEXT*)(exts +k);      /* and occurrence extensions */
  for (i = 0; i < k; i++) {     /* initialize the pattern extensions */
    e = exts+i; e->supp = 0; e->cnt = 0; e->oxs = x; x += rd.frqs[i]; }
  for (j = 0; j < n; j++) {     /* traverse the transactions and */
    o = occs +j;                /* the items in each transaction */
    for (s = o->items; s->item >= 0; s++) {
      e = exts   +s->item;      /* get the corresp. pattern extension */
      x = e->oxs +e->cnt++;     /* and its next occurrence extension */
      x->item  = s;             /* set the extension item */
      x->occ   = o;             /* and the pattern occurrence */
      e->supp += o->wgt;        /* sum transaction weights (support) */
    }                           /* (exts represents the possible */
  }                             /* extensions of the empty sequence) */
  memset(rd.frqs, 0, (size_t)k *sizeof(TID));
  r = rec_iw(exts, z, 0, &rd);  /* search for frequent sequences */
  if ((r < tbg_wgt(tabag))      /* report empty sequence if closed */
  || !(mode & ISR_CLOSED))      /* or all sequences are requested */
    r = (isr_isetx(report, NULL, 0, NULL, tbg_wgt(tabag), 0, 0) < 0)
      ? -1 : 0;                 /* report the empty sequence */
  free(rd.wgts);                /* deallocate the pattern buffer, */
  free(exts); free(occs);       /* pattern extensions and occurrences */
  return (r < 0) ? (int)r : 0;  /* return the error status */
}  /* sequoia_iw() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef SOIA_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("item weight output format characters (option -i#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%w  sum of the item weights\n");
  printf("  %%m  mean/average of the item weights\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%i  number of items (sequence length)\n");
  printf("  %%a  absolute sequence support\n");
  printf("  %%s  relative sequence support as a fraction\n");
  printf("  %%S  relative sequence support as a percentage\n");
  printf("  %%Q  total transaction weight (database size)\n");
  printf("All format characters can be preceded by the number\n");
  printf("of significant digits to be printed (at most 32 digits),\n");
  printf("even though this value is ignored for integer numbers.\n");
  #endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

/*--------------------------------------------------------------------*/

static ITEM getbdr (char *s, char **end, double **border)
{                               /* --- get the support border */
  ITEM   i, k;                  /* loop variables */
  double *b;                    /* support border */

  assert(s && end && border);   /* check the function arguments */
  for (i = k = 0; s[i]; i++)    /* traverse the string and */
    if (s[i] == ':') k++;       /* count the number separators */
  *border = b = (double*)malloc((size_t)++k *sizeof(double));
  if (!b) return -1;            /* allocate a support border */
  for (i = 0; i < k; i++) {     /* traverse the parameters */
    b[i] = strtod(s, end);      /* get the next parameter and */
    if (*end == s) break;       /* check for an empty parameter */
    s = *end; if (*s++ != ':') break;
  }                             /* check for a colon */
  if (++i < k)                  /* shrink support array if possible */
    *border = (double*)realloc(b, (size_t)i *sizeof(double));
  return i;                     /* return number of support values */
}  /* getbdr() */

/*--------------------------------------------------------------------*/

static int setbdr (ISREPORT *report, SUPP w, ITEM zmin,
                   double **border, ITEM n)
{                               /* --- set the support border */
  double s;                     /* to traverse the support values */

  assert(report                 /* check the function arguments */
  &&    (w > 0) && (zmin >= 0) && border && (*border || (n <= 0)));
  while (--n >= 0) {            /* traverse the support values */
    s = (*border)[n];           /* transform to absolute count */
    s = ceilsupp((s >= 0) ? s/100.0 *(double)w *(1-DBL_EPSILON) : -s);
    if (isr_setbdr(report, n+zmin, (RSUPP)s) < 0) return -1;
  }                             /* set support in item set reporter */
  if (*border) { free(*border); *border = NULL; }
  return 0;                     /* return 'ok' */
}  /* setbdr() */

/*--------------------------------------------------------------------*/

#ifndef NDEBUG                  /* if debug version */
  #undef  CLEANUP               /* clean up memory and close files */
  #define CLEANUP \
  if (twrite) twr_delete(twrite, 1); \
  if (report) isr_delete(report, 0); \
  if (tabag)  tbg_delete(tabag,  0); \
  if (tread)  trd_delete(tread,  1); \
  if (ibase)  ib_delete (ibase);     \
  if (border) free(border);
#endif

GENERROR(error, exit)           /* generic error reporting function */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- main function */
  int     i, k = 0;             /* loop variables, counters */
  char    *s;                   /* to traverse the options */
  CCHAR   **optarg = NULL;      /* option argument */
  CCHAR   *fn_inp  = NULL;      /* name of input  file */
  CCHAR   *fn_out  = NULL;      /* name of output file */
  CCHAR   *fn_psp  = NULL;      /* name of pattern spectrum file */
  CCHAR   *recseps = NULL;      /* record  separators */
  CCHAR   *fldseps = NULL;      /* field   separators */
  CCHAR   *wgtseps = NULL;      /* weight  separators */
  CCHAR   *blanks  = NULL;      /* blank   characters */
  CCHAR   *comment = NULL;      /* comment characters */
  CCHAR   *hdr     = "";        /* record header  for output */
  CCHAR   *sep     = " ";       /* item separator for output */
  CCHAR   *iwf     = ":%m";     /* output format for item weights */
  CCHAR   *dflt    = " (%S)";   /* default format for check */
  CCHAR   *info    = dflt;      /* format for information output */
  int     target   = 's';       /* target type (closed/maximal) */
  ITEM    zmin     = 1;         /* minimum length of a sequence */
  ITEM    zmax     = ITEM_MAX;  /* maximum length of a sequence */
  double  supp     = 10;        /* minimum support (in percent) */
  SUPP    smin     = 1;         /* minimum support of an item set */
  int     mtar     = TA_DUPERR; /* mode for transaction reading */
  int     scan     = 0;         /* flag for scanable item output */
  int     bdrcnt   = 0;         /* number of support values in border */
  int     stats    = 0;         /* flag for sequence statistics */
  PATSPEC *psp;                 /* collected pattern spectrum */
  ITEM    m;                    /* number of items */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  #ifndef QUIET                 /* if not quiet version */
  clock_t t;                    /* timer for measurements */

  prgname = argv[0];            /* get program name for error msgs. */

  /* --- print usage message --- */
  if (argc > 1) {               /* if arguments are given */
    fprintf(stderr, "%s - %s\n", argv[0], DESCRIPTION);
    fprintf(stderr, VERSION); } /* print a startup message */
  else {                        /* if no arguments given */
    printf("usage: %s [options] infile [outfile]\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-t#      target type                              "
                    "(default: %c)\n", target);
    printf("         (s: frequent, c: closed sequences)\n");
    printf("-m#      minimum number of items per sequence     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per sequence     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of a sequence            "
                    "(default: %g%%)\n", supp);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-F#:#..  support border for filtering item sets   "
                    "(default: none)\n");
    printf("         (list of minimum support values, "
                    "one per item set size,\n");
    printf("         starting at the minimum size, "
                    "as given with option -m#)\n");
    printf("-P#      write a pattern spectrum to a file\n");
    printf("-Z       print item set statistics "
                    "(number of item sets per size)\n");
    printf("-g       write output in scanable form "
                    "(quote certain characters)\n");
    printf("-h#      record header  for output                "
                    "(default: \"%s\")\n", hdr);
    printf("-k#      item separator for output                "
                    "(default: \"%s\")\n", sep);
    printf("-i#      output format for item weights           "
                    "(default: \"%s\")\n", iwf);
    printf("-v#      output format for sequence information   "
                    "(default: \"%s\")\n", info);
    printf("-w       integer transaction weight in last field "
                    "(default: only items)\n");
    printf("-r#      record/transaction separators            "
                    "(default: \"\\n\")\n");
    printf("-f#      field /item        separators            "
                    "(default: \" \\t,\")\n");
    printf("-b#      blank   characters                       "
                    "(default: \" \\t\\r\")\n");
    printf("-u#      weight  separators                       "
                    "(default: none)\n");
    printf("-C#      comment characters                       "
                    "(default: \"#\")\n");
    printf("-!       print additional option information\n");
    printf("infile   file to read transactions from           "
                    "[required]\n");
    printf("outfile  file to write frequent sequences to      "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: acdejlopqwxyz [A-Z]\[CFPZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 's';      break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': supp   =       strtod(s, &s);    break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'g': scan   = 1;                      break;
          case 'h': optarg = &hdr;                   break;
          case 'k': optarg = &sep;                   break;
          case 'i': optarg = &iwf;                   break;
          case 'v': optarg = &info;                  break;
          case 'w': mtar  |= TA_WEIGHT;              break;
          case 'r': optarg = &recseps;               break;
          case 'f': optarg = &fldseps;               break;
          case 'b': optarg = &blanks;                break;
          case 'u': optarg = &wgtseps;               break;
          case 'C': optarg = &comment;               break;
          default : error(E_OPTION, *--s);           break;
        }                       /* set option variables */
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       /* get option argument */
    else {                      /* -- if argument is no option */
      switch (k++) {            /* evaluate non-options */
        case  0: fn_inp = s;      break;
        case  1: fn_out = s;      break;
        default: error(E_ARGCNT); break;
      }                         /* note filenames */
    }
  }
  if (optarg)       error(E_OPTARG);     /* check option arguments */
  if (k      < 1)   error(E_ARGCNT);     /* and number of arguments */
  if (zmin   < 0)   error(E_SIZE, zmin); /* check the size limits */
  if (zmax   < 0)   error(E_SIZE, zmax); /* and the minimum support */
  if (supp   > 100) error(E_SUPPORT, supp);
  if (bdrcnt < 0)   error(E_NOMEM);
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'c': target = ISR_CLOSED;           break;
    /* case 'm': target = ISR_MAXIMAL;          break; */
    default : error(E_TARGET, (char)target); break;
  }                             /* (get the target type code) */
  if (info == dflt)             /* adapt the default info. format */
    info = (supp < 0) ? " (%a)" : " (%S)";
  MSG(stderr, "\n");            /* terminate the startup message */

  /* --- read transaction database --- */
  ibase = ib_create((wgtseps && *wgtseps) ? IB_WEIGHTS : 0, 0);
  if (!ibase) error(E_NOMEM);   /* create an item base */
  tabag = tbg_create(ibase);    /* and a transaction bag */
  if (!tabag) error(E_NOMEM);   /* to store the transactions */
  tread = trd_create();         /* create a table reader and */
  if (!tread) error(E_NOMEM);   /* configure the characters */
  trd_allchs(tread, recseps, fldseps, blanks, "", comment);
  if (wgtseps && *wgtseps) {    /* if there are weight separators */
    trd_chars(tread, TA_WGTSEP,          wgtseps);
    trd_chars(tread, TRD_FLDSEP|TRD_ADD, wgtseps);
  }                             /* set them as add. field separators */
  CLOCK(t);                     /* start timer, open input file */
  if (trd_open(tread, NULL, fn_inp) != 0)
    error(E_FOPEN, trd_name(tread));
  MSG(stderr, "reading %s ... ", trd_name(tread));
  k = tbg_read(tabag, tread, mtar);
  if (k < 0) error(-k, tbg_errmsg(tabag, NULL, 0));
  trd_delete(tread, 1);         /* read the transaction database, */
  tread = NULL;                 /* then delete the table scanner */
  m = ib_cnt(ibase);            /* get the number of items, */
  n = tbg_cnt(tabag);           /* the number of transactions, */
  w = tbg_wgt(tabag);           /* the total transaction weight */
  MSG(stderr, "[%"ITEM_FMT" item(s), %"TID_FMT, m, n);
  if (w != (SUPP)n) { MSG(stderr, "/%"SUPP_FMT, w); }
  MSG(stderr, " transaction(s)] done [%.2fs].", SEC_SINCE(t));
  if ((m <= 0) || (n <= 0))     /* check for at least one item */
    error(E_NOITEMS);           /* and at least one transaction */
  MSG(stderr, "\n");            /* compute absolute support value */
  supp = (supp >= 0) ? supp/100.0 *(double)w *(1-DBL_EPSILON) : -supp;
  smin = (SUPP)ceilsupp(supp);  /* compute absolute support values */

  /* --- sort and recode items --- */
  CLOCK(t);                     /* start timer, print log message */
  MSG(stderr, "recoding items ... ");
  m = tbg_recode(tabag, 0, -1, -1, -1);
  if (m <  0) error(E_NOMEM);   /* recode items and transactions */
  if (m <= 0) error(E_NOITEMS); /* check the number of items */
  MSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- trim and reduce transactions --- */
  CLOCK(t);                     /* start timer, print log message */
  MSG(stderr, "filtering and reducing transactions ... ");
  tbg_sort(tabag, 1, 0);        /* sort the trans. lexicographically */
  n = tbg_reduce(tabag, 0);     /* and reduce them to unique ones */
  MSG(stderr, "[%"TID_FMT, n);  /* print number of transactions */
  if (w != (SUPP)n) { MSG(stderr, "/%"SUPP_FMT, w); }
  MSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));

  /* --- find frequent sequences --- */
  report = isr_create(ibase);   /* create an item set reporter */
  if (!report) error(E_NOMEM);  /* and configure it */
  isr_setsize(report, zmin, zmax);
  isr_setsupp(report, smin, SUPP_MAX);
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmtx(report, scan, hdr, sep, NULL, info, iwf) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report)); /* open the item set file */
  if ((isr_settarg(report, target, ISR_NOFILTER, -1) < 0)
  ||  (isr_setup(report) < 0))  /* set target and reporting mode */
    error(E_NOMEM);             /* and set up the item set reporter */
  CLOCK(t);                     /* start timer, print log message */
  MSG(stderr, "writing %s ... ", isr_name(report));
  k = (wgtseps && *wgtseps)     /* search for frequent sequences */
    ? sequoia_iw(tabag, target, smin, 0, report)
    : sequoia   (tabag, target, smin, 0, report);
  if (k < 0) error(E_NOMEM);    /* check for a search error */
  MSG(stderr, "[%"SIZE_FMT" sequence(s)]", isr_repcnt(report));
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  if (stats)                    /* print item set statistics */
    isr_prstats(report, stdout, 0);
  if (isr_close(report) != 0)   /* close the output file */
    error(E_FWRITE, isr_name(report));

  /* --- write pattern spectrum --- */
  if (fn_psp) {                 /* if to write a pattern spectrum */
    CLOCK(t);                   /* start timer, create table write */
    psp    = isr_getpsp(report);/* get the pattern spectrum */
    twrite = twr_create();      /* create a table writer and */
    if (!twrite) error(E_NOMEM);/* open the output file */
    if (twr_open(twrite, NULL, fn_psp) != 0)
      error(E_FOPEN,  twr_name(twrite));
    MSG(stderr, "writing %s ... ", twr_name(twrite));
    if (psp_report(psp, twrite, 1.0) != 0)
      error(E_FWRITE, twr_name(twrite));
    twr_delete(twrite, 1);      /* write the pattern spectrum */
    twrite = NULL;              /* and delete the table writer */
    MSG(stderr, "[%"SIZE_FMT" signature(s)]", psp_sigcnt(psp));
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* write a log message */

  /* --- clean up --- */
  CLEANUP;                      /* clean up memory and close files */
  SHOWMEM;                      /* show (final) memory usage */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif
