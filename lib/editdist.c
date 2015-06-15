/* 
Contains 3 primary functions

### PRIMER ALIGNMENT
bounded_edit_distance - edit distance bounded to the 5' end with restrictions to the 3' end

edit_distance - standard edit distance

### BARCODE ALIGNMENT
hammingdist_distance - hamming distance 
*/

#include "Python.h"

#if defined(_MSC_VER)
typedef unsigned __int8 u_int8_t;
#endif

#define EDITDIST_VERSION    "0.6"

/* $Id: editdist.c,v 1.5 2007/05/03 23:36:36 djm Exp $ */
/* $Id: editdist.c,v 0.4 2013/12/31 mls $ */
/* $Id: editdist.c,v 0.4 2014/4/17 mls added hamming distnace */
/* $Id: editdist.c,v 0.6 2015/6/13 mls change to allow for faster processing */

#ifndef MIN
# define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

typedef struct Tuple {
    int dist;
    int pos;
} Tuple;

/* 
compute the Levenstein distance between a (primer) and b (sequence read)
pegged to the 5' end and bounded by edit distance k
*/
struct Tuple
bounded_edit_distance(const char *a, int alen, const char *b, int blen, int k, int m)
{
    // a is primer, b is seq, k is max error and m is end matches
    int i, j;
    int *current, *previous, *tmpl, add, del, chg;
    Tuple val = { k+1, alen };
    /* a (primer) should always be < b (read) */
    if (alen > blen) {
        return (val);
    }

    if (alen == 0){
        val.dist = -2;
        return (val);
    }
    if ((previous = calloc(alen + 1, sizeof(*previous))) == NULL) {
        free(previous);
        val.dist = -1;
        return (val);
    }
    if ((current = calloc(alen + 1, sizeof(*current))) == NULL) {
        free(current);
        val.dist = -1;
        return (val);
    }

    for (i = 0; i < alen - m + 1; i++)
        previous[i] = i;

    for (i = 1; i < alen - m + k + 1 ; i++) {
        if (i > 1) {
            memset(previous, 0, (alen - m + 1) * sizeof(*previous));
            tmpl = previous;
            previous = current;
            current = tmpl;
        }
        current[0] = i;
        for (j = 1; j < alen -m + 1; j++) {
            add = previous[j] + 1;
            del = current[j - 1] + 1;
            chg = previous[j - 1];
            if (a[j - 1] != b[i - 1])
                chg++;
            current[j] = MIN(add, del);
            current[j] = MIN(current[j], chg);
        }
        if (current[alen - m] <= val.dist ){
            val.dist = current[alen - m ]; // bottom right node, global alignment
            val.pos = i+m;
        }
    }

    for (i = 1; i <= m; i++){
        if (a[alen - i] != b[val.pos - i]){
            val.dist = val.dist+100; // penelty of 100 for not meeting end match criteria
            break;          
        }       
    }
    free(previous);
    free(current);

    return (val);
}


PyDoc_STRVAR(bounded_editdist_distance_doc,
"bounded_edit_distance(a, b, k, m) -> int, int\n\
    Calculates the bounded Levenshtein's edit distance between strings \"a\" and \"b\" with bound \"k\" and \"m\" matching bases at end\n");

static PyObject *
bounded_editdist_distance(PyObject *self, PyObject *args)
{
    Tuple r;
    char *a, *b;
    int alen, blen, k, m;

    if (!PyArg_ParseTuple(args, "s#s#ii", &a, &alen, &b, &blen, &k, &m))
                return NULL;
    r = bounded_edit_distance(a, alen, b, blen, k, m);
    if (r.dist== -1) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return NULL;
    }
    if (r.dist== -2) {
        PyErr_SetString(PyExc_SystemError, "Bad Arguments");
        return NULL;
    }
    return Py_BuildValue("ii", r.dist, r.pos);
}

/*
Compute the Levenstein distance between a and b 
*/
static int
edit_distance(const char *a, int alen, const char *b, int blen)
{
    int tmplen, i, j;
    const char *tmp;
    int *current, *previous, *tmpl, add, del, chg, r;

    /* Swap to reduce worst-case memory requirement */
    if (alen > blen) {
        tmp = a;
        a = b;
        b = tmp;
        tmplen = alen;
        alen = blen;
        blen = tmplen;
    }

    if (alen == 0)
        return (blen);

    if ((previous = calloc(alen + 1, sizeof(*previous))) == NULL) {
        free(previous);
        return (-1);
    }
    if ((current = calloc(alen + 1, sizeof(*current))) == NULL) {
        free(current);
        return (-1);
    }
    for (i = 0; i < alen + 1; i++)
        previous[i] = i;

    for (i = 1; i < blen + 1; i++) {
        if (i > 1) {
            memset(previous, 0, (alen + 1) * sizeof(*previous));
            tmpl = previous;
            previous = current;
            current = tmpl;
        }
        current[0] = i;
        for (j = 1; j < alen + 1; j++) {
            add = previous[j] + 1;
            del = current[j - 1] + 1;
            chg = previous[j - 1];
            if (a[j - 1] != b[i - 1])
                chg++;
            current[j] = MIN(add, del);
            current[j] = MIN(current[j], chg);
        }
    }
    r = current[alen];
    free(previous);
    free(current);
    return (r);
}

PyDoc_STRVAR(editdist_distance_doc,
"distance(a, b) -> int\n\
    Calculates Levenshtein's edit distance between strings \"a\" and \"b\"\n");

static PyObject *
editdist_distance(PyObject *self, PyObject *args)
{
    char *a, *b;
    int alen, blen, r;

    if (!PyArg_ParseTuple(args, "s#s#", &a, &alen, &b, &blen))
                return NULL;
    r = edit_distance(a, alen, b, blen);
    if (r == -1) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return NULL;
    }
    return PyInt_FromLong(r);
}

/*
Compute the hamming distance between a and b 
*/
static int
hammingdist_distance(const char *a, int alen, const char *b, int blen, int score)
{
    int i, diff;

    /* a and b should be same length, otherwise return error */
    if (alen != blen){
        return (-2);
    }
    /* length should not be 0, otherwise return error */
    if (alen == 0)
        return (-2);

    diff = 0;
    for (i = 0; i < alen; i++){
        if (a[i] != b[i])
            diff++;
        if (diff > score){
            return (diff);
        }
    }
    return (diff);
}

PyDoc_STRVAR(hamming_distance_doc,
"distance(a, b) -> int\n\
    Calculates hamming distance between two equal length strings \"a\" and \"b\"\n");

static PyObject *
hamming_distance(PyObject *self, PyObject *args)
{
    char *a, *b;
    int alen, blen, r;

    if (!PyArg_ParseTuple(args, "s#s#", &a, &alen, &b, &blen))
                return NULL;

    r = hammingdist_distance(a, alen, b, blen, blen);
    if (r == -2) {
        PyErr_SetString(PyExc_SystemError, "Bad Arguments");
        return NULL;
    }
    return PyInt_FromLong(r);
}

PyDoc_STRVAR(hamming_distance_list_doc,
"distance(a, b) -> int\n\
    Calculates hamming distance between a string and a list of strings\n");

static PyObject *
hamming_distance_list(PyObject *self, PyObject *args)
{
    char *b;
    int blen, b, s, r;
    PyListObject* barcode_list = NULL;

    if (!PyArg_ParseTuple(args, "O!s#", &PyList_Type, &barcode_list_o, &b, &blen)){
        return NULL;
    }

    b = 0; // set the initial best barcode to the first
    s = blen; // set the initial score to the length of the barcode

    Py_ssize_t barcode_list_o_length = PyList_GET_SIZE(barcode_list_o);
    for (Py_ssize_t i = 0; i < barcode_list_o_length; i++) {
        PyObject * barcode_o = PyList_GET_ITEM(barcode_list_o, i);
        if (PyString_Check(barcode_o)) {
            r = hammingdist_distance(PyString_AS_STRING(barcode_o), (int)PyString_GET_SIZE(barcode_o), b, blen, s);
            if (r == -2) {
                PyErr_SetString(PyExc_SystemError, "Bad Arguments");
                return NULL;
            } else if ( r < s ){
                b = (int)i;
                s = r;
            }
        } else {
            PyErr_SetString(PyExc_TypeError,
                "first argument must be a list of strings");
            return NULL;
        }
    }
    Py_DECREF(barcode_list);
    return Py_BuildValue("ii", b, s);
}

static PyMethodDef editdist_methods[] = {
    {   "distance", (PyCFunction)editdist_distance,
        METH_VARARGS,   editdist_distance_doc       },
    {   "bounded_distance", (PyCFunction)bounded_editdist_distance,
        METH_VARARGS,   bounded_editdist_distance_doc       },
    {   "hamming_distance", (PyCFunction)hamming_distance,
        METH_VARARGS,    hamming_distance_doc       },
    {   "hamming_distance_list", (PyCFunction)hamming_distance_list,
        METH_VARARGS,    hamming_distance_list_doc  },
    {   NULL, NULL, 0, NULL }  /* sentinel */
};

PyDoc_STRVAR(module_doc, "Calculate Hamming distance, Levenshtein's edit distance, and a edge bounded Levenshtein's edit distance.\n");

PyMODINIT_FUNC
initeditdist(void)
{
    PyObject *m;

    m = Py_InitModule3("editdist", editdist_methods, module_doc);
    PyModule_AddStringConstant(m, "__version__", EDITDIST_VERSION);
}
