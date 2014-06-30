#! /usr/bin/env python

from __future__ import print_function

import re
import os


re_single_line_comment = re.compile( r'//.*$', flags = re.MULTILINE )
re_multi_line_comment = re.compile( r'/[*].*?[*]/', flags = re.DOTALL )
re_multi_whitespace = re.compile( r'\s+' )
re_whitespace_comma = re.compile( r'\s*,\s?' )
re_whitespace_closing_parenthesis = re.compile( r'\s*\)\s*' )
re_whitespace_opening_parenthesis = re.compile( r'\s*\(\s*' )
re_export = re.compile( r'EXPORT .*?\(.*?\)', flags = re.DOTALL )


def parse( filename, funcdefs ):

    with open( filename, 'r' ) as file:
        file = file.read()
    # remove all comments
    file, n = re_single_line_comment.subn( lambda m: '', file )
    file, n = re_multi_line_comment.subn( lambda m: '', file )
    # find EXPORT and write to `funcdefs`
    for match in re_export.findall( file ):
        match = match[7:] + ';'
        match, n = re_multi_whitespace.subn( lambda m: ' ', match )
        match, n = re_whitespace_comma.subn( lambda m: ', ', match )
        match, n = re_whitespace_closing_parenthesis.subn( lambda m: ' )', match )
        match, n = re_whitespace_opening_parenthesis.subn( lambda m: '( ', match )
        print( match, file = funcdefs )


if __name__ == '__main__':

    import sys

    source_files = []
    output = None

    args = iter( sys.argv )
    next( args ) # skip program name
    for arg in args:
        if arg == '--':
            break
        elif arg == '--output':
            output = next( args )
        else:
            source_files.append( arg )
    source_files.extend( args )

    if output is None:
        print( 'ERROR: `--output OUTPUT` argument missing' )
        raise SystemExit( 1 )

    with open( output, 'w' ) as funcdefs:
        print( '#pragma once', file = funcdefs )
        print( '#define EXPORT', file = funcdefs )
        path_to_headers = os.path.relpath( 'src/headers', os.path.dirname( output ) )
        print( '#include "{0}/header_clsses.h"'.format( path_to_headers ), file = funcdefs )
        print( '#include "{0}/array.h"'.format( path_to_headers ), file = funcdefs )
        for file in source_files:
            if not file.endswith( '.cpp' ):
                continue
            parse( file, funcdefs )
