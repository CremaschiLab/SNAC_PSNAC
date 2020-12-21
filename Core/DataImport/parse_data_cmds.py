
#  _________________________________________________________________________
#
#  Coopr: A COmmon Optimization Python Repository
#  Copyright (c) 2008 Sandia Corporation.
#  This software is distributed under the BSD License.
#  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
#  the U.S. Government retains certain rights in this software.
#  For more information, see the Coopr README.txt file.
#  _________________________________________________________________________

__all__ = ['parse_data_commands']

import sys
import re
import os
import os.path
import ply.lex as lex
import ply.yacc as yacc
from pyutilib.misc import flatten_list
from pyutilib.ply import t_newline, t_ignore, _find_column, p_error, ply_init

## -----------------------------------------------------------
##
## Lexer definitions for tokenizing the input
##
## -----------------------------------------------------------

_parse_info = None
debugging = False

reserved = {
    'data' : 'DATA',
    'set' : 'SET',
    'param' : 'PARAM',
    'end' : 'END',
    'import' : 'IMPORT',
    'include' : 'INCLUDE',
    'namespace' : 'NAMESPACE',
    }

# Token names
tokens = [
    "COMMA",
    "LBRACE",
    "RBRACE",
#    "NUMBER",
    "SEMICOLON",
    "COLON",
    "COLONEQ",
    "LBRACKET",
    "RBRACKET",
    "LPAREN",
    "RPAREN",
#    "RANGE",
    "WORD",
    "WORDWITHINDEX",
    "WORDWITHSQUOTEDINDEX",
    "STRING",
    "QUOTEDSTRING",
    "FILENAME",
    "EQ",
    "TR",
    "ASTERISK",
    "NONWORD",
] + list(reserved.values())

# Regular expression rules
t_COMMA     = r","
t_LBRACKET  = r"\["
t_RBRACKET  = r"\]"
t_LBRACE  = r"\{"
t_RBRACE  = r"\}"
#t_NUMBER    = r"[0-9]+(\.[0-9]+){0,1}"
t_SEMICOLON = r";"
t_COLON     = r":"
t_COLONEQ   = r":="
t_EQ        = r"="
t_TR        = r"\(tr\)"
#t_LT        = r"<"
#t_GT        = r">"
t_LPAREN    = r"\("
t_RPAREN    = r"\)"
t_ASTERISK  = r"\*"

# Discard comments
def t_COMMENT(t):
    r'\#[^\n]*'
    #global _comment_list
    #_comment_list.append(t.value)

def t_WORDWITHINDEX(t):
    r'[a-zA-Z_0-9][a-zA-Z_0-9\.\-]*\[[a-zA-Z_0-9\.\-,\*]*\]'
    if t.value in reserved:
        t.type = reserved[t.value]    # Check for reserved words
    #t.type = reserved.get(t.value,'WORDWITHINDEX')    # Check for reserved words
    return t

def t_WORDWITHSQUOTEDINDEX(t):
    r'[a-zA-Z_0-9][a-zA-Z_0-9\.\-]*\[[\'\"a-zA-Z_0-9\.\-,\* ]*\]'
    if t.value in reserved:
        t.type = reserved[t.value]    # Check for reserved words
    #t.type = reserved.get(t.value,'WORDWITHSQUOTEDINDEX')    # Check for reserved words
    return t

def t_WORD(t):
    r'[a-zA-Z_0-9][a-zA-Z_0-9\.+\-]*'
    if t.value in reserved:
        t.type = reserved[t.value]    # Check for reserved words
    #t.type = reserved.get(t.value,'WORD')    # Check for reserved words
    return t

def t_STRING(t):
    r'[a-zA-Z_0-9\.+\-]+'
    if t.value in reserved:
        t.type = reserved[t.value]    # Check for reserved words
    #t.type = reserved.get(t.value,'STRING')    # Check for reserved words
    return t

def t_QUOTEDSTRING(t):
    r'"([^"]|\"\")*"|\'([^\']|\'\')*\''
    if t.value in reserved:
        t.type = reserved[t.value]    # Check for reserved words
    #t.type = reserved.get(t.value,'QUOTEDSTRING')    # Check for reserved words
    return t

def t_FILENAME(t):
    r'[a-zA-Z_0-9\./\\]*(/|\\)[a-zA-Z_0-9\-\./\\]*|[a-zA-Z_0-9\./\\]*(/|\\)[a-zA-Z_0-9\-\./\\]*'
    if t.value in reserved:
        t.type = reserved.get(t.value)    # Check for reserved words
    else:
        t.type = 'FILENAME'
    #t.type = reserved.get(t.value,'FILENAME')    # Check for reserved words
    return t

t_NONWORD   = r"[^\.A-Za-z0-9,;:=<>\*\(\)\#{}\[\] \n\t\r]+"

# Error handling rule
def t_error(t):             #pragma:nocover
    raise IOError
    t.lexer.skip(1)


## -----------------------------------------------------------
##
## Yacc grammar for data commands
##
## -----------------------------------------------------------

def p_expr(p):
    '''expr : statements
            | '''
    if len(p) == 2:
        #print "STMTS",p[1]
        for stmt in p[1]:
            if type(stmt) is list:
                _parse_info[None].append(stmt)
            else:
                for key in stmt:
                    if key in _parse_info:
                        _parse_info[key].append(stmt[key])
                    else:
                        _parse_info[key] = stmt[key]

def p_statements(p):
    '''statements : statements statement
                  | statement
                  | statements NAMESPACE WORD LBRACE statements RBRACE
                  | NAMESPACE WORD LBRACE statements RBRACE '''
    #print "STMT X",p[1:],p[1]
    len_p = len(p)
    if len_p == 3:
        # NB: statements will never be None, but statement *could* be None
        p[0] = p[1]
        if p[2] is not None:
            p[0].append(p[2])
    elif len_p == 2:
        if p[1] is None:
            p[0] = []
        else:
            p[0] = [p[1]]
    elif len_p == 7:
        # NB: statements will never be None
        p[0] = p[1]
        p[0].append({p[3]:p[5]})
    else:
        # NB: statements will never be None
        p[0] = [{p[2] : p[4]}]

def p_statement(p):
    '''statement : SET WORD COLONEQ setdecl SEMICOLON
                 | SET WORD COLONEQ SEMICOLON
                 | SET WORD COLON items COLONEQ setdecl SEMICOLON
                 | SET WORD COLON items COLONEQ SEMICOLON
                 | SET WORDWITHINDEX COLONEQ setdecl SEMICOLON
                 | SET WORDWITHINDEX COLONEQ SEMICOLON
                 | SET WORDWITHSQUOTEDINDEX COLONEQ setdecl SEMICOLON
                 | SET WORDWITHSQUOTEDINDEX COLONEQ SEMICOLON
                 | PARAM items COLONEQ paramdecl SEMICOLON
                 | IMPORT importdecl SEMICOLON
                 | INCLUDE WORD SEMICOLON
                 | INCLUDE QUOTEDSTRING SEMICOLON
                 | DATA SEMICOLON
                 | END SEMICOLON
    '''
    #print "STATEMENT",len(p), p[1:]
    stmt = p[1]
    if stmt == 'set' or stmt == 'param':
        p[0] = flatten_list(p[1:-1])
    elif stmt == 'include':
        p[0] = p[1:-1]
    elif stmt == 'import':
        p[0] = [p[1]]+ p[2]
    else:
        # Not necessary, but nice to document how statement could end up None
        p[0] = None 

def p_setdecl(p):
    '''setdecl : items'''
    p[0] = p[1]

def p_paramdecl(p):
    '''paramdecl : items'''
    p[0] = p[1]

def p_importdecl(p):
    '''importdecl : filename import_options
                  | filename
                  | filename import_options COLON WORD EQ indices variable_options
                  | filename COLON WORD EQ indices variable_options
                  | filename import_options COLON indices variable_options
                  | filename COLON indices variable_options
                  | filename import_options COLON variable_options
                  | filename COLON variable_options
    '''
    tmp = {'filename':p[1]}
    if len(p) == 2:
        p[0] = [tmp, (None,[]), {}]
    elif len(p) == 3:
        tmp.update(p[2])
        p[0] = [tmp, (None,[]), {}]
    elif len(p) == 4:
        p[0] = [tmp, (None,[]), p[3]]
    elif len(p) == 5:
        if p[2] == ':':
            p[0] = [tmp, (None,p[3]), p[4]]
        else:
            tmp.update(p[2])
            p[0] = [tmp, (None,[]), p[4]]
    elif len(p) == 6:
        tmp.update(p[2])
        p[0] = [tmp, (None,p[4]), p[5]]
    elif len(p) == 7:
        p[0] = [tmp, (p[3],p[5]), p[6]]
    elif len(p) == 8:
        tmp.update(p[2])
        p[0] = [tmp, (p[4],p[6]), p[7]]
    else:
        raise IOError

def p_import_options(p):
    '''import_options : WORD EQ STRING import_options
                      | WORD EQ STRING
                      | WORD EQ QUOTEDSTRING import_options
                      | WORD EQ QUOTEDSTRING
                      | WORD EQ WORD import_options
                      | WORD EQ WORD
                      | WORD EQ PARAM import_options
                      | WORD EQ PARAM
                      | WORD EQ SET import_options
                      | WORD EQ SET
    '''
    tmp = {p[1]:p[3]}
    if len(p) == 4:
        p[0] = tmp
    else:
        tmp.update(p[4])
        p[0] = tmp

def p_variable_options(p):
    '''variable_options : variable variable_options
                        | variable
    '''
    if len(p) == 2:
        p[0] = p[1]
    else:
        p[1].update(p[2])
        p[0] = p[1]

def p_variable(p):
    '''variable : WORD
                | WORD EQ WORD
    '''
    if len(p) == 2:
        p[0] = {p[1]:p[1]}
    else:
        p[0] = {p[3]:p[1]}

def p_indices(p):
    '''indices : LBRACKET WORD index_list RBRACKET
               | LBRACKET WORD RBRACKET
    '''
    if len(p) == 5:
        p[0] = p[3]
        p[0].insert(0,p[2])
    else:
        p[0] = [p[2]]

def p_index_list(p):
    '''index_list : COMMA WORD index_list
                  | COMMA ASTERISK index_list
                  | COMMA WORD
                  | COMMA ASTERISK
    '''
    if len(p) == 4:
        p[0] = p[3]
        p[0].insert(0,p[2])
    else:
        p[0] = [p[2]]

def p_set_template(p):
    '''set_template : LPAREN WORD index_list RPAREN
                | LPAREN ASTERISK index_list RPAREN
                | LPAREN WORD RPAREN
                | LPAREN ASTERISK RPAREN
    '''
    if len(p) == 5:
        p[0] = p[1]+",".join([p[2]]+p[3])+p[4]
    else:
        p[0] = p[1]+p[2]+p[3]

def p_param_template(p):
    '''param_template : LBRACKET WORD index_list RBRACKET
                | LBRACKET ASTERISK index_list RBRACKET
                | LBRACKET WORD RBRACKET
                | LBRACKET ASTERISK RBRACKET
    '''
    if len(p) == 5:
        p[0] = p[1]+",".join([p[2]]+p[3])+p[4]
    else:
        p[0] = p[1]+p[2]+p[3]

# Composite p_item into p_items
# Produces ~3% speedup, but at the cost of readability
def p_items(p):
    '''
    items : items WORD
          | items WORDWITHINDEX
          | items WORDWITHSQUOTEDINDEX
          | items NONWORD
          | items STRING
          | items QUOTEDSTRING
          | items COMMA
          | items COLON
          | items LBRACE
          | items RBRACE
          | items LBRACKET
          | items RBRACKET
          | items TR
          | items LPAREN
          | items RPAREN
          | items ASTERISK
          | items set_template
          | items param_template
          | WORD
          | WORDWITHINDEX
          | WORDWITHSQUOTEDINDEX
          | NONWORD
          | STRING
          | QUOTEDSTRING
          | COMMA
          | COLON
          | LBRACE
          | RBRACE
          | LBRACKET
          | RBRACKET
          | TR
          | LPAREN
          | RPAREN
          | ASTERISK
          | set_template
          | param_template
    '''
    # Locate and handle item as necessary
    single_item = len(p) == 2
    if single_item:
        tmp = p[1]
    else:
        tmp = p[2]
    if tmp[0] == '"' and tmp[-1] == '"' and len(tmp) > 2 and not ' ' in tmp:
        tmp = tmp[1:-1]

    # Grow items list according to parsed item length
    if single_item:
        p[0] = [tmp]
    else:
        # yacc __getitem__ is expensive: use a local list to avoid a
        # getitem call on p[0]
        tmp_lst = p[1]
        tmp_lst.append(tmp)
        p[0] = tmp_lst

def p_filename(p):
    '''filename : WORD
                | STRING
                | QUOTEDSTRING
                | FILENAME
                | WORD COLON FILENAME
    '''
    if len(p) == 2:
        p[0] = p[1]
    else:
        p[0] = p[1]+p[2]+p[3]

#
# the ampl dat file lexer and yaccer only need to be
# created once, so have the corresponding objects
# accessible at module scope.
#

tabmodule = 'parse_table_datacmds'

ampl_dat_lexer = None
ampl_dat_yaccer = None

#
# The function that performs the parsing
#
def parse_data_commands(data=None, filename=None, debug=0, outputdir=None):

    global debugging
    global ampl_dat_lexer
    global ampl_dat_yaccer

    if outputdir is None:
        outputdir = os.path.dirname(os.path.abspath(__file__))+os.sep

    # if the lexer/yaccer haven't been initialized, do so.
    if ampl_dat_lexer is None:
        #
        # Always remove the parser.out file, which is generated to create debugging
        #
        if os.path.exists("parser.out"):        #pragma:nocover
            os.remove("parser.out")
        if debug > 0:                           #pragma:nocover
            #
            # Remove the parsetab.py* files.  These apparently need to be removed
            # to ensure the creation of a parser.out file.
            #
            if os.path.exists("parsetab.py"):
                os.remove("parsetab.py")
            if os.path.exists("parsetab.pyc"):
                os.remove("parsetab.pyc")
            debugging=True

        ampl_dat_lexer = lex.lex()
        #
        tmpsyspath = sys.path
        sys.path.append(outputdir)
        ampl_dat_yaccer = yacc.yacc(debug=debug, tabmodule=tabmodule, outputdir=outputdir)
        sys.path = tmpsyspath

    #
    # Initialize parse object
    #
    global _parse_info
    _parse_info = {}
    _parse_info[None] = []

    #
    # Parse the file
    #
    global _parsedata
    if not data is None:
        _parsedata=data
        ply_init(_parsedata)
        ampl_dat_yaccer.parse(data, lexer=ampl_dat_lexer, debug=debug)
    elif not filename is None:
        f = open(filename, 'r')
        try:
            data = f.read()
        except Exception:
            f.close()
            del f
            raise e
        f.close()
        del f
        _parsedata=data
        ply_init(_parsedata)
        ampl_dat_yaccer.parse(data, lexer=ampl_dat_lexer, debug=debug)
    else:
        _parse_info = None
    #
    # Disable parsing I/O
    #
    debugging=False
    return _parse_info
