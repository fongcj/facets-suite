#!/opt/common/CentOS_6-dev/python/python-2.7.10/bin/python

import argparse, os, sys, re
import subprocess, itertools, errno, csv, gzip
import imp, glob
## import cmo

FACETS_DIR = None


def make_sure_path_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def gzip_file_with_size(file_path):
    if not os.path.exists(file_path): return(False)
    with gzip.open(file_path, 'rb') as f:
        for i, l in enumerate(f):
            i = i+1
            if i > 10: return(True)
        return(False)


def slugify(value):
    """
    Normalizes string, removes non-alpha characters,
    and converts spaces to hyphens.
    http://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename-in-python
    """
    ##    import unicodedata
    ##    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore')
    value = re.sub('[^\w\s-]', '', value).strip()
    re.sub('[-\s]+', '-', value)
    return(value)



def run(args, facets_args):
    try:
        import cmo
        os.environ['PYTHON2']=cmo.util.programs['python']['default']
    except:
        os.environ['PYTHON2']="/usr/bin/env python"
    cmd = [args.script]
    del vars(args)['script']
    for key, value in vars(args).items():
        if value == True:
            cmd = cmd + ['--'+key]
        elif value == False:
            pass
            #don't append
        elif value != None:
            if type(value) is list:
                cmd = cmd + ["--"+key] +  value
            else:
                cmd = cmd + ["--"+key, str(value)]
    print >>sys.stderr, "Executing %s" % " ".join(cmd)
    try:
        rv = subprocess.call(cmd + facets_args)
        if rv!=0:
            sys.exit(rv)
    except:
        # hack: when facets are running part of Roslin pipeline,
        # any kind of facets failure prevents downstream workflows from running
        # due to how toil works.
        # this workaround reports exitcode 0 to toil even if there is an error
        if os.environ.get('FACETS_OVERRIDE_EXITCODE'):
            sys.exit(0)
        else:
            print >>sys.stderr, "Error executing command, returncode %d" % rv
            sys.exit(rv)


def add_subparser(file, subparsers):
    R_script = open(file, "r")
    arguments_found = 0
    finished = 1
    parser = None
    for line in R_script:
        #begin perversion of all that is right and holy
        #take parser$add_argument lines out of R scripts
        #strip out the "type" field, which doesn't match pythons datatypes
        #and interpret the result as a line of python code
        #we do this in the name of good, i.e. so our subparsers can reproduce the help
        #but what have we really wrought?
        #behold I am become monkey, mixer of languages
        add_argument_index = line.find("add_argument")
        if finished== 0:
            python_code += line.rstrip()
            if line.find(")") > -1:
                finished = 1
                #boosh R to python
                exec python_code
        if add_argument_index > -1:
            arguments_found = 1
            if not parser:
                #parser name will be filename without .R
                parser = subparsers.add_parser(name=os.path.basename(file).replace(".R",""))
            #take "add_argument" to end of line
            python_code = line[add_argument_index:].rstrip()
            #take out "type" argument from R
            fixed_code =re.sub("""type\s*=\s*['"]\S+['"],""","", python_code)
            #probably not needed unless someone is a jackass
            fixed_code = re.sub("""required\s*=\s*T""", """required=True""", fixed_code)
            fixed_code = re.sub("""required\s*=\s*F""", """required=False""", fixed_code)
            fixed_code = re.sub("""default\s*=\s*F(alse|ALSE)?\s*,""", """""", fixed_code)
            fixed_code = re.sub("""default\s*=\s*T(rue|RUE)?\s*,""", """""", fixed_code)
            python_code =  "parser." +fixed_code
            if python_code.find(")") == -1:
                #multiline shit
                finished=0
            else:
                #boosh R to python
                exec python_code
    if parser:
        parser.set_defaults(script=file)
    #extra special sprinkle of shit-magic to get subparser help to print if no args
    if len(sys.argv) ==2 and sys.argv[1] == os.path.basename(file).replace(".R",""):
        parser.print_help()
        sys.exit(1)



FACETS_DIR = os.path.dirname(os.path.realpath(__file__))
def create_parser(parser=None):
    if not parser:
        parser = argparse.ArgumentParser(description="run FACETS analysis")
    subparsers = parser.add_subparsers(help='sub-command help')
    for file in glob.glob(os.path.join(FACETS_DIR, "*.R")):
        add_subparser(file, subparsers)
    return parser




if __name__ =='__main__':

    ### ARGUMENTS
    parser = create_parser()
    (args, facets_args) = parser.parse_known_args()
    args_dict= vars(args)
    run(args, facets_args)
##    print args


