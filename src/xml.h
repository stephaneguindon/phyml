#include <config.h>

#ifndef XML_H
#define XML_H

#include "utilities.h"

xml_node *XML_Load_File(FILE *fp);
int XML_Add_Character(int c, char  **bufptr, char **buffer, int *bufsize);
int XML_Parse_Element(FILE *fp, xml_node *n);
int XML_Set_Attribute(xml_node *n, char *attr_name, char *attr_value);
xml_attr *XML_Make_Attribute(xml_attr *prev, char *attr_name, char *attr_value);
void XML_Make_Node_Id(xml_node *n, char *id);
int XML_Set_Node_Id(xml_node *n, char *id);
int XML_Set_Node_Value(xml_node *n, char *val);
void XML_Make_Node_Value(xml_node *n, char *val);
xml_node *XML_Search_Node_Name(char *name, int skip, xml_node *node);
char *XML_Get_Attribute_Value(xml_node *node, char *id);
int XML_Validate_Attr_Int(char *target, int num, ...);
void XML_Free_XML_Attr(xml_attr *attr);
void XML_Free_XML_Node(xml_node *node);
void XML_Free_XML_Tree(xml_node *node);
xml_node *XML_Search_Node_ID(char *id, int skip, xml_node *node);
xml_node *XML_Make_Node(char *name);
xml_node *XML_Search_Node_Attribute_Value(char *attr_name, char *value, int skip, xml_node *node);
void XML_Check_Siterates_Node(xml_node *prev);
int XML_Get_Number_Of_Classes_Siterates(xml_node *prev);
int XML_Siterates_Number_Of_Classes(xml_node *sr_node);
void XML_Check_Duplicate_ID(xml_node *n);
void XML_Count_Number_Of_Node_With_ID(char *id, int *count, xml_node *n);
void XML_Count_Number_Of_Node_With_Name(char *name, int *count, xml_node *n);
void XML_Write_XML_Graph(FILE *fp, xml_node *root);
void XML_Write_XML_Node(FILE *fp, int *indent, xml_node *root);
xml_attr *XML_Search_Attribute(xml_node *n, char *target_attr_name);
xml_node *XML_Copy_XML_Graph(xml_node *root);
void XML_Copy_XML_Node(xml_node *cpy_root, xml_node *root);
void Check_Mandatory_XML_Node(xml_node *root, char *name);
xml_node *XML_Search_Node_Generic(char *nd_name, char *attr_name, char *attr_val, int skip, xml_node *node);


#endif
