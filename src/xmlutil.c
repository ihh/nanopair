#include <string.h>
#include "xmlutil.h"
#include "util.h"

/* private builder methods */
xmlNode* getNodeByName (xmlNode* node, const char* name) {
  xmlNode* child;
  Assert (node != NULL, "XML parse error while searching for <%s>\n", name);
  for (child = node->children; child; child = child->next)
    if (child->type == XML_ELEMENT_NODE && strcmp ((const char*) child->name, name) == 0)
      return child;
  Abort ("XML parse error while searching in <%s>: <%s> not found\n", node->name, name);
  return (xmlNode*) NULL;
}

xmlChar* getNodeContent (xmlNode* node) {
  Assert (node != NULL, "XML parse error\n");
  return (xmlChar*) node->children->content;
}

xmlChar* getNodeContentOrComplain (xmlNode* node, const char* tag) {
  if (node == NULL) {
    if (tag)
      Abort ("Node is null, expected <%s>\n", tag);
    else
      Abort ("Node and tag are null\n");
  }
  if (node->children == NULL)
    Abort ("Missing children for tag: %s\n", tag);
  return (xmlChar*) node->children->content;
}

xmlChar* getAttrByName (xmlNode* node, const char* name) {
  xmlAttr* attr;
  Assert (node != NULL, "XML parse error while searching for attribute %s\n", name);
  for (attr = node->properties; attr; attr = attr->next)
    if (strcmp ((const char*) attr->name, name) == 0)
      return (xmlChar*) attr->content;
  return (xmlChar*) NULL;
}
