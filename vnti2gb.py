#!/usr/bin/env python
# encoding: utf-8
'''
 vnti2gb -- Command line tool for converting VectorNTI *.mol and *.seq files into GenBank records

@author:     Aaron Kitzmiller

@copyright:  2017 The Presidents and Fellows of Harvard College. All rights reserved.

@license:    GPL v2.0

@contact:   aaron_kitzmiller@harvard.edu
@deffield    updated: Updated
'''

import sys, os, traceback, logging, re

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

DEBUG = os.environ.get('VNTI2GB_DEBUG', False)
# Setup logging
loglevel = logging.INFO
if DEBUG:
    loglevel = logging.DEBUG

logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s", level=loglevel)
logger = logging.getLogger('vnti2gb')

# Setup regex matches
regexes = {
    'LOCUS'         : re.compile(r'^25\|(.*)'),
    'LENGTH'        : re.compile(r'^33\|(.*)'),
    'DEFINITION'    : re.compile(r'^34\|(.*)'),
    'SOURCE'        : re.compile(r'^1007\|(.*)'),
    'DATE'          : re.compile(r'^1006\|(.*)'),
    'REFERENCE'     : re.compile(r'^1009\|(.*)'),
    'ORGANISM1'     : re.compile(r'^1015\|(.*)'),
    'ORGANISM2'     : re.compile(r'^1018\|(.*)'),
}
# Regexes for features
fregexes = {
    'LABEL'         : re.compile(r'^52\|(.*)'),
    'LOCSTR'        : re.compile(r'^285\|(.*)'),
    'FKEY'          : re.compile(r'^51\|(.*)'),
    'NOTE'          : re.compile(r'^54\|(.*)'),
    'TAGS'          : re.compile(r'^286\|(.*)'),
    'START'         : re.compile(r'^55\|(.*)'),
    'END'           : re.compile(r'^56\|(.*)'),
}
# Feature keys
fkeys = {
    '4': 'CDS',
    '5': 'centromere',
    '9': 'enhancer',
    '12': 'homeodomain',
    '13': 'iDNA',
    '14': 'insertion_seq',
    '15': 'intron',
    '18': 'loci',
    '19': 'repeat_region',
    '20': 'misc_binding',
    '21': 'misc_feature',
    '22': 'misc_marker',
    '23': 'modified_base',
    '25': 'regulatory',
    '26': 'polyA_site',
    '27': 'primer',
    '28': 'primer_bind',
    '29': 'promoter',
    '30': 'promoter',
    '31': 'protein_bind',
    '32': 'RBS',
    '33': 'rep_origin',
    '34': 'repeat_region',
    '35': 'repeat_region',
    '38': 'splicing_signal',
    '40': 'STS',
    '41': 'TATA_signal',
    '43': 'terminator',
    '44': 'transposon',
    '46': 'ZF_domain',
    '47': '-10_signal',
    '48': '-35_signal',
    '50': "3'UTR",
    '51': "5'clip",
    '52': "5'UTR",
    '53': 'misc_RNA',
    '54': 'mRNA',
    '56': 'scRNA',
    '58': 'rRNA',
    '59': 'tRNA',
    '60': 'gene',
    '61': 'exon',
    '62': 'exon',
    '81': 'allele',
    '83': 'misc_difference',
    '84': 'mat_peptide',
    '85': 'misc_difference',
    '86': 'misc_recomb',
    '87': 'misc_signal',
    '88': 'misc_structure',
    '89': 'old_sequence',
    '95': 'transit_peptide',
    '94': 'sig_peptide',
    '96': 'snp',
    '97': 'virion',
    '98': 'source',
    '99': 'unsure',
    '102': 'gap',
    '274': 'overhang',
}

# Reverse the keys so that I can assign /vntifkey
fkeyrev = dict((v,k) for k,v in fkeys.iteritems())

# Pad for beginning of feature tags
tagpadding = 21 * ' '


def makeGenBankReport(molfile,seqfile,commentfile=None):
    rpt = '''{LOCUS}
{DEFINITION}
{SOURCE}
{REFERENCE}
{COMMENT}
{FEATURES}
{SEQUENCE}
'''
    data = {}
    data['DEFINITION'] = ''
    data['SOURCE'] = ''
    data['ORGANISM1'] = ''

    # Get the sequence
    with open(seqfile,'r') as f:
        seq = ''.join(f.readlines())
    if seq.strip() == '':
        raise UserException('No sequence data in %s' % seqfile)
    data['SEQ'] = seq
    data['FEATURES'] = []

    # Go through the molfile a line at a time
    lines = []
    with open(molfile,'r') as f:
        lines = f.readlines()
    if len(lines) == 0:
        raise UserException('No lines in mol file %s' % molfile)

    # Get the comment lines
    commentlines = []
    if commentfile:
        with open(commentfile,'r') as f:
            for line in f:
                parts = line.split('|')
                if len(parts) > 1:
                    commentlines.append(parts[1].strip())
        if len(commentlines) == 0:
            raise UserException('No lines in comment file %s' % commentfile)



    featuredata = None
    for line in lines:
        # End of a feature
        if line.strip() == '50' and featuredata is not None:
            # End of a feature
            data['FEATURES'].append(featuredata)
            featuredata = None
        if line.strip() == '45':
            # Starting a new feature
            featuredata = {}

        if featuredata is not None:
            # check the feature regexes
            for field, regex in fregexes.iteritems():
                m = regex.match(line)
                if m:
                    if field == 'FKEY':
                        # If it's a feature key, do the translation from number to name
                        if m.group(1).strip() in fkeys:
                            featuredata[field] = fkeys[m.group(1).strip()]
                        else:
                            logger.info('Unable to match feature key %s' % str(m.group(1)))
                    else:
                        featuredata[field] = m.group(1).strip()

                    continue

        # Check single line regexes
        for field, regex in regexes.iteritems():
            m = regex.match(line)
            if m:
                data[field] = m.group(1).strip()
                continue


    # Compose the report
    # LOCUS line
    locusstart = '{LABEL:12}{LOCUS}'.format(LABEL='LOCUS',LOCUS=data['LOCUS'])
    locusextra = '{LENGTH} bp    DNA    circular    {DATE}'.format(LENGTH=data['LENGTH'],DATE='01-JAN-2017')
    pad = 72 - (len(locusstart) + len(locusextra))
    if pad < 1:
        pad = 1
    LOCUS = locusstart + pad * ' ' + locusextra

    # DEFINITION line
    DEFINITION = '{LABEL:12}{DEFINITION:.70}'.format(LABEL='DEFINITION',DEFINITION=data['DEFINITION'])

    # SOURCE line(s)
    SOURCE = '{LABEL:12}{SOURCE:.70}\n  {SUBLABEL}  {ORGANISM1:.70}'.format(
        LABEL='SOURCE',SOURCE=data['SOURCE'],SUBLABEL='ORGANISM',ORGANISM1=data['ORGANISM1']
    )
    if 'ORGANISM2' in data:
        SOURCE += '\n' + 12 * ' ' + data['ORGANISM2']

    # REFERENCE section
    REFERENCE = ''
    if 'REFERENCE' in data:
        REFERENCE = data['REFERENCE'].replace('|','\n')

    # COMMENT section
    COMMENT = ''
    if commentlines:
        commentjoin = '\n' + ' ' * 12
        COMMENT = '{LABEL:12}{COMMENTS}\n'.format(LABEL='COMMENT',COMMENTS=commentjoin.join(commentlines))

    # Features
    FEATURES = ''
    if len(data['FEATURES']) > 0:
        FEATURES = 'FEATURES' + 13 * ' ' + 'Location/Qualifiers\n'
    for feature in data['FEATURES']:
        if 'FKEY' not in feature:
            logger.info('Unable to create feature due to missing key')
        elif 'LOCSTR' not in feature and 'START' not in feature and 'END' not in feature:
            logger.info('Unable to create feature due to missing location information')
        else:
            # Setup location string
            locstr = ''
            if 'LOCSTR' in feature:
                locstr = feature['LOCSTR']
            elif 'START' in feature and 'END' in feature:
                locstr = '{START}..{END}'.format(START=feature['START'],END=feature['END'])
            else:
                raise Exception('Not sure how you got here')

            featurestr = '     {FKEY:16}{LOCSTR}\n'.format(FKEY=feature['FKEY'],LOCSTR=locstr)
            if 'TAGS' in feature:
                tags = feature['TAGS'].split('|')
                tagjoin = '\n' + tagpadding
                featurestr += tagpadding + tagjoin.join(tags) + '\n'
            if 'LABEL' in feature:
                featurestr += tagpadding + '/label="%s"\n' % feature['LABEL']
            if 'NOTE' in feature:
                featurestr += tagpadding + '/note="%s"\n' % feature['NOTE']
            featurestr += tagpadding + '/vntifkey="%s"\n' % fkeyrev[feature['FKEY']]
            FEATURES += featurestr


    # The sequence
    # Get counts
    seqstr = data['SEQ'].lower()
    counts = dict((b,seqstr.count(b)) for b in ['a','c','g','t'])
    SEQUENCE = 'BASE COUNT     %d a      %d c      %d g      %d t\n' % (counts['a'], counts['c'], counts['g'], counts['t'])
    SEQUENCE += 'ORIGIN\n'

    # Split into 60 char lines
    seqlines = [seqstr[i:i + 60] for i in range(0, len(seqstr), 60)]
    for i, line in enumerate(seqlines):
        # split into 10 char chunks
        chunks = [line[j:j + 10] for j in range(0, 60, 10)]
        SEQUENCE += '{POS:>9} {SEQ}\n'.format(POS=str(i * 60 + 1), SEQ=' '.join(chunks))

    

    return rpt.format(
        LOCUS=LOCUS,
        DEFINITION=DEFINITION,
        SOURCE=SOURCE,
        REFERENCE=REFERENCE,
        COMMENT=COMMENT,
        FEATURES=FEATURES,
        SEQUENCE=SEQUENCE,
    )


class UserException(Exception):
    '''
    I can actually get a message from this exception
    '''
    def __init__(self, message):
        super(UserException, self).__init__(message)
        self.user_msg = message


def main(argv=None):

    program_name = 'vnti2gb'
    program_version_message = '0.1'
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by Aaron Kitzmiller on %s.
  Copyright 2017 The Presidents and Fellows of Harvard College. All rights reserved.

  Licensed under GPL v2.0
  http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, '1/9/2017')

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('dbdir',metavar="DBDIR",help="Directory containing the VectorNTI database")
        parser.add_argument('outputdir',metavar='OUTPUTDIR',help='Genbank output directory. Files will be named to match up with *.mol and *.seq')

        # Process arguments
        args = parser.parse_args()

        vntidbpath = args.dbdir
        if not os.path.exists(vntidbpath):
            raise UserException('Path, %s, does not exist' % vntidbpath)

        outputpath = args.outputdir
        if not os.path.exists(outputpath):
            raise UserException('Output dir, %s, does not exist' % outputpath)


        # Read file names from MolData, match up Seq file and process
        moldatapath = os.path.join(vntidbpath,'MolData')
        seqdatapath = os.path.join(moldatapath,'Seq')
        commentdatapath = os.path.join(moldatapath,'Comment')
        molfiles = os.listdir(moldatapath)
        seqfiles = os.listdir(seqdatapath)
        commentfiles = os.listdir(commentdatapath)

        logger.info('%d files in MolData dir %s' % (len(molfiles),moldatapath))
        logger.info('%d files in Seq dir %s' % (len(seqfiles),seqdatapath))
        logger.info('%d files in Comment dir %s' % (len(commentfiles),commentdatapath))

        for molfile in molfiles:
            if os.path.isfile(os.path.join(moldatapath,molfile)):
                logger.debug('Processing mol file %s' % molfile)
                seqfile = molfile.replace('mol','seq')
                commentfile = molfile.replace('mol','cmn')
                if seqfile in seqfiles:
                    try:
                        if commentfile in commentfiles:
                            gbstr = makeGenBankReport(
                                os.path.join(moldatapath,molfile),
                                os.path.join(seqdatapath,seqfile),
                                os.path.join(commentdatapath,commentfile),
                            )
                        else:
                            gbstr = makeGenBankReport(
                                os.path.join(moldatapath,molfile),
                                os.path.join(seqdatapath,seqfile)
                            )
                        gbpath = os.path.join(outputpath,molfile.replace('mol','gb'))
                        with open(gbpath,'w') as gb:
                            gb.write(gbstr)
                        logger.debug('Wrote %s' % gbpath)
                    except Exception as e:
                        logger.info('Unable to make report for mol file %s: %s' % (molfile,str(e)))
                else:
                    logger.info('Unable to find seq file %s, to match mol file %s' % (seqfile,molfile))


    except KeyboardInterrupt:
        # handle keyboard interrupt #
        return 0
    except Exception as e:
        if hasattr(e, 'user_msg') and not DEBUG:
            sys.stderr.write("%s: %s\n" % (program_name, e.user_msg))
        else:
            sys.stderr.write("%s: %s\n%s\n" % (program_name, str(e),traceback.format_exc()))
        return 2


if __name__ == "__main__":
    sys.exit(main())
