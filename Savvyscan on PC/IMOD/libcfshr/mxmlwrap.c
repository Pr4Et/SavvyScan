/*
 * imodxml.c - An index-based interface to mxml library
 *
 * Copyright (C) 2016 by the Regents of the University of 
 * Colorado.  See dist/COPYRIGHT for full notice.
 *
 * $Id$
 */                                                                           

#include "b3dutil.h"
#include "ilist.h"
#include "parse_params.h"
#include "mxmlwrap.h"
#include "mxml.h"

#define NODE_LIST_QUANTUM  32

static Ilist **sNodeLists = NULL;
static int sNumLists = 0;
static int sLastLevel = -1;

static int getOrAddFreeList();
static mxml_node_t *getNodeAtIndex(int xmlInd, int nodeInd, int *error);
static int getElementString(int xmlInd, int nodeInd, const char **string);

/*!
 * Reads an XML file from [filename] and returns a copy of the name of its root element in
 * [rootElement], which should be freed with 
 * {free}.  The node index list is initialized with a pointer to the top xml node 
 * at index 0 and a pointer to the root element at index 1.  The return value is the 
 * index for this set of XML data, or -1 for a memory error, -2 for failure to open file,
 * or -3 for failure to load as XML.
 */
int ixmlReadFile(const char *filename, char **rootElement)
{
  FILE *fp;
  int xmlInd;
  mxml_node_t *xml, *top;
  const char *root;
  
  fp = fopen(filename, "r");
  if (!fp)
    return -2;
  xml = mxmlLoadFile(NULL, fp, MXML_OPAQUE_CALLBACK);
  if (!xml)
    return -3;
  fclose(fp);

  xmlInd = getOrAddFreeList();
  if (xmlInd < 0) {
    mxmlDelete(xml);
    return xmlInd;
  }
  ilistAppend(sNodeLists[xmlInd], &xml);

  /* Get the top node and then get the next if it is opaque node */
  top = mxmlWalkNext(xml, xml, MXML_DESCEND);
  if (top->type == MXML_OPAQUE)
    top = mxmlWalkNext(top, xml, MXML_DESCEND);
  ilistAppend(sNodeLists[xmlInd], &top);
  root = mxmlGetElement(top);
  *rootElement = NULL;
  if (root)
    *rootElement = strdup(root);
  return xmlInd;
}

/*!
 * Starts a new set of XML data and a new node list, with a root element named by
 * [rootElement].  The node index list is initialized with a pointer to the top xml node 
 * at index 0 and a pointer to the root element at index 1.  The return value is the 
 * index for this set of XML data, -1 for a memory error, or -2 
 * for failure to create a node.
 */
int ixmlNewNodeList(const char *rootElement)
{
  mxml_node_t *xml, *top;
  int xmlInd = getOrAddFreeList();
  if (xmlInd >= 0) {
    xml = mxmlNewXML("1.0");
    if (!xml)
      xmlInd = -2;
    else {
      ilistAppend(sNodeLists[xmlInd], &xml);
      top = mxmlNewElement(xml, rootElement);
      if (!top)
        xmlInd = -2;
      else
        ilistAppend(sNodeLists[xmlInd], &top);
    }
  }
  return xmlInd;
}

/*!
 * Writes the XML data at index [xmlInd] to file with name [filename].
 */
int ixmlWriteFile(int xmlInd, const char *filename)
{
  FILE *afile;
  int err;
  mxml_node_t *xml = getNodeAtIndex(xmlInd, 0, &err);
  if (!xml)
    return err;
  afile = fopen(filename, "w");
  if (!afile)
    return -1;
  mxmlSetWrapMargin(0);
  sLastLevel = -1;
  err = mxmlSaveFile(xml, afile, (mxml_save_cb_t)ixmlWhitespace_cb);
  fclose(afile);
  return err;
}

/*!
 * Looks for all elements with the name [tag] under the node with list index [nodeInd] in
 * the XML data at index [xmlInd], and adds pointers to the elements to the node list.
 * The number of elements found is returned in [numFound] and the list index of the first 
 * one found is returned in [foundInd].  The return value is -1 for a memory error,
 * -2 for an invalid XML index, or -3 for an invalid node index.
 */
int ixmlFindElements(int xmlInd, int nodeInd, const char *tag, int *foundInd, 
                     int *numFound)
{
  mxml_node_t *node, *parent;
  const char *key;
  int err;
  *numFound = 0;
  parent = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!parent)
    return err;
  node = mxmlGetFirstChild(parent);
  *foundInd = ilistSize(sNodeLists[xmlInd]);
  while (node != NULL) {
    key = mxmlGetElement(node);
    if (key && !strcmp(key, tag)) {
      if (ilistAppend(sNodeLists[xmlInd], &node))
        return -1;
      (*numFound)++;
    }
    node = mxmlGetNextSibling(node);
  }
  return 0;
}

/*!
 * Returns a copy of the element value for the node with list index [nodeInd] in
 * the XML data at index [xmlInd] into [string], which should be freed with {free}.
 * The return value is -1 for a memory error, -2 for an invalid XML index, -3 for an 
 * invalid node index, -4 for an invalid value node, such as multiple child nodes or
 * wrong kind of node.
 */
int ixmlGetStringValue(int xmlInd, int nodeInd, char **string)
{
  const char *value;
  int err = getElementString(xmlInd, nodeInd, &value);
  if (err)
    return err;
  *string = strdup(value);
  if (!*string)
    return -1;
  return 0;
}

/*!
 * Returns an integer value for the node with list index [nodeInd] in the XML data at 
 * index [xmlInd] into [val].  The return value is -1 for errors in parsing the value,
 * -2 for an invalid XML index, -3 for an invalid node index, or -4 for an invalid value 
 * node.
 */
int ixmlGetIntegerValue(int xmlInd, int nodeInd, int *val)
{
  const char *string;
  int numToGet = 1;
  int err = getElementString(xmlInd, nodeInd, &string);
  if (err)
    return err;
  return PipGetLineOfValues(string, string, (void *)val, PIP_INTEGER, &numToGet, 1);
}

/*!
 * Like @@ixmlGetIntegerValue@, but returns a floating point value.
 */
int ixmlGetFloatValue(int xmlInd, int nodeInd, float *val)
{
  const char *string;
  int numToGet = 1;
  int err = getElementString(xmlInd, nodeInd, &string);
  if (err)
    return err;
  return PipGetLineOfValues(string, string, (void *)val, PIP_FLOAT, &numToGet, 1);
}

/*!
 * Like @@ixmlGetIntegerValue@, but returns a double value.
 */
int ixmlGetDoubleValue(int xmlInd, int nodeInd, double *val)
{
  const char *string;
  int numToGet = 1;
  int err = getElementString(xmlInd, nodeInd, &string);
  if (err)
    return err;
  return PipGetLineOfValues(string, string, (void *)val, PIP_DOUBLE, &numToGet, 1);
}

/*!
 * Returns a copy of the attribute with the name [name] for the node with list index 
 * [nodeInd] in the XML data at index [xmlInd] into [string], which should be freed with
 * {free}.  The return value is -1 for a memory error, -2 for an invalid XML index,
 * or -3 for an invalid node index.
 */
int ixmlGetStringAttribute(int xmlInd, int nodeInd, const char *name, char **string)
{
  int err, ind;
  mxml_node_t *node;
  node = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!node)
    return err;
  for (ind = 0; ind < node->value.element.num_attrs; ind++) {
    if (!strcmp(node->value.element.attrs[ind].name, name)) {
      *string = strdup(node->value.element.attrs[ind].value);
      if (!*string)
        return -1;
      return 0;
    }
  }
  return 1;
}

/*!
 * Returns the integer value of the attribute with the name [name] for the node with list 
 * index [nodeInd] in the XML data at index [xmlInd] into [val].
 * The return value is -1 for errors in parsing the value,
 * -2 for an invalid XML index, or -3 for an invalid node index.
 */
int ixmlGetIntegerAttribute(int xmlInd, int nodeInd, const char *name, int *val)
{
  int err, ind;
  mxml_node_t *node;
  int numToGet = 1;
  node = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!node)
    return err;
  for (ind = 0; ind < node->value.element.num_attrs; ind++)
    if (!strcmp(node->value.element.attrs[ind].name, name))
      return PipGetLineOfValues(node->value.element.attrs[ind].value, 
                                node->value.element.attrs[ind].value,
                                (void *)val, PIP_INTEGER, &numToGet, 1);
  return 1;
}

/*!
 * Adds an element with the name [tag] to the node with list index [nodeInd] in the XML 
 * data at index [xmlInd].  The new element is added to the node index list.  Returns the
 * list index  of the node, or -1 for a memory error, -2 for an invalid 
 * XML index, or -3 for an invalid node index.
 */
int ixmlAddElement(int xmlInd, int nodeInd, const char *tag)
{
  int err;
  mxml_node_t *node, *elem;
  node = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!node)
    return err;
  elem = mxmlNewElement(node, tag);
  if (!elem)
    return -1;
  if (ilistAppend(sNodeLists[xmlInd], &elem))
    return -1;
  return ilistSize(sNodeLists[xmlInd]) - 1;
}

/*!
 * Sets the value of the node with list index [nodeInd] in the XML
 * data at index [xmlInd] to [string].  Returns -1 for a memory error, -2 for an invalid
 * XML index, or -3 for an invalid node index.
 */
int ixmlSetStringValue(int xmlInd, int nodeInd, const char *string)
{
  int err = 0;
  mxml_node_t *node;
  node = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!node)
    return err;
  if (!mxmlNewText(node, 0, string))
    err = -1;
  return err;
}

/*!
 * Like @@ixmlSetStringValue@, but sets the integer value in [val].
 */
int ixmlSetIntegerValue(int xmlInd, int nodeInd, int val)
{
  char buffer[64];
  sprintf(buffer, "%d", val);
  return ixmlSetStringValue(xmlInd, nodeInd, buffer);
}

/*!
 * Like @@ixmlSetStringValue@, but sets the float value in [val].
 */
int ixmlSetFloatValue(int xmlInd, int nodeInd, float val)
{
  char buffer[64];
  sprintf(buffer, "%g", val);
  return ixmlSetStringValue(xmlInd, nodeInd, buffer);
}

/*!
 * Like @@ixmlSetStringValue@, but sets the double value in [val].
 */
int ixmlSetDoubleValue(int xmlInd, int nodeInd, double val)
{
  char buffer[64];
  sprintf(buffer, "%g", val);
  return ixmlSetStringValue(xmlInd, nodeInd, buffer);
}

/*!
 * Adds a string attribute with the given [name] and [value] to the node with list index 
 * [nodeInd] in the XML data at index [xmlInd].
 * Returns -2 for an invalid XML index, or -3 for an invalid node index.
 */
int ixmlAddStringAttribute(int xmlInd, int nodeInd, const char *name, const char *value)
{
  int err;
  mxml_node_t *node;
  node = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!node)
    return err;
  mxmlElementSetAttr(node, name, value);
  return 0;
}

/*!
 * Like @@ixmlAddStringAttribute@, but sets the value string from the integer in [value].
 */
int ixmlAddIntegerAttribute(int xmlInd, int nodeInd, const char *name, int value)
{
  int err;
  mxml_node_t *node;
  node = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!node)
    return err;
  mxmlElementSetAttrf(node, name, "%d", value);
  return 0;
}

/*!
 * Discards the end of the node list starting at index [firstInd] in the in the XML data 
 * at index [xmlInd].
 */
int ixmlPopListIndexes(int xmlInd, int firstInd)
{
  if (xmlInd < 0 || xmlInd >= sNumLists || !sNodeLists[xmlInd])
    return -2;
  if (firstInd < 2 || firstInd >= ilistSize(sNodeLists[xmlInd]))
    return -3;
  ilistTruncate(sNodeLists[xmlInd], firstInd);
  return 0;
}

/*!
 * Clears out the XML structure and index list for XML item [index].
 */
void ixmlClear(int index)
{
  mxml_node_t **xml;
  if (index < 0 && index >= sNumLists)
    return;
  xml = ilistItem(sNodeLists[index], 0);
  mxmlDelete(*xml);
  ilistDelete(sNodeLists[index]);
  sNodeLists[index] = NULL;
}

/*!
 * Resets the output level counter before writing to file with @ixmlWhitespace_cb
 */
void ixmlResetLastLevel()
{
  sLastLevel = -1;
}

/*!
 * Callback function for writing with 2 spaces of indentation, must be case to
 * mxml_save_cb_t in the call to mxmlSaveFile.
 * Returns a Whitespace string or NULL; [nodeVoid] is an Element node, [where] is a
 * MXML_WS_ tag.
 */
const char *ixmlWhitespace_cb(void *nodeVoid, int where)
{
  mxml_node_t	*parent;		/* Parent node */
  int		level;			/* Indentation level */
  /* Tabs for indentation */
  /* static const char *tabs = "\t\t\t\t\t\t\t\t"; */
  char spaces[33] = {"                                "};
  static char buffer[36];
  mxml_node_t *node = (mxml_node_t *)nodeVoid;

  if (where != MXML_WS_BEFORE_OPEN && where != MXML_WS_BEFORE_CLOSE)
    return NULL;
  for (level = -1, parent = node->parent;
       parent;
       level ++, parent = parent->parent);
  if (level > 16)
    level = 16;
  else if (level < 0)
    level = 0;
  if (sLastLevel < 0) {
    sLastLevel = level;
    return NULL;
  }
  if (level == sLastLevel && where == MXML_WS_BEFORE_CLOSE)
    return NULL;
  sLastLevel = level;
  sprintf(buffer, "\n%s", &spaces[32 - 2 * level]);
  return &buffer[0];
}

/*
 * Creates a new Ilist structure and places it pointer in a free spot of the array
 * or extends the array, returns the index of the new list or -1 for memory errors.
 */
static int getOrAddFreeList()
{
  Ilist **newLists;
  int xmlInd;
  for (xmlInd = 0; xmlInd < sNumLists; xmlInd++)
    if (!sNodeLists[xmlInd])
      break;
  if (xmlInd >= sNumLists) {
    newLists = (Ilist **)malloc((sNumLists + 1) * sizeof(Ilist *));
    if (!newLists)
      return -1;
    if (sNumLists)
      memcpy(newLists, sNodeLists, sNumLists * sizeof(Ilist *));
    free(sNodeLists);
    sNodeLists = newLists;
    sNumLists++;
  }
  sNodeLists[xmlInd] = ilistNew(sizeof(mxml_node_t *), NODE_LIST_QUANTUM);
  if (!sNodeLists[xmlInd])
    xmlInd = -1;
  else
    ilistQuantum(sNodeLists[xmlInd], NODE_LIST_QUANTUM);
  return xmlInd;
} 

/*
 * Makes sure xmlInd and nodeInd are legal and if so, returns the node pointer at those
 * indexes.  Otherwise returns NULL and sets error to -2 or -3.
 */
static mxml_node_t *getNodeAtIndex(int xmlInd, int nodeInd, int *error)
{
  Ilist *nodeList;
  mxml_node_t **listPtr;
  if (xmlInd < 0 || xmlInd >= sNumLists) {
    *error = -2;
    return NULL;
  }
  nodeList = sNodeLists[xmlInd];
  if (!nodeList) {
    *error = -2;
    return NULL;
  }
  if (nodeInd < 0 || nodeInd >= ilistSize(nodeList)) {
    *error = -3;
    return NULL;
  }
  *error = 0;
  listPtr = (mxml_node_t **)ilistItem(nodeList, nodeInd);
  return *listPtr;
}

/*
 * Gets the original pointer an the element value string at the given indexes.
 * Returns the usual -2 or -3 errors, or -4 if there is not one child node or it
 * is not an opaque node.
 */
static int getElementString(int xmlInd, int nodeInd, const char **string)
{
  int err;
  mxml_node_t *node, *child;
  node = getNodeAtIndex(xmlInd, nodeInd, &err);
  if (!node)
    return err;
  child = mxmlGetFirstChild(node);

  if (!child || mxmlGetType(child) != MXML_OPAQUE || mxmlGetLastChild(node) != child)
    return -4;
  *string = child->value.opaque;
  return 0;
}
