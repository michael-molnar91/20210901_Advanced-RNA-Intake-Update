from functools import reduce
from model_utils import Choices
from buckaneer.product.exceptions import SequenceValidationException
from buckaneer.product.utils import get_sequence_length_range
from buckaneer.synthego.common.oligolib.chemistry.sequence import Sequence as ChemistrySequence
from buckaneer.product.constants import PRODUCT_SLUG_CHOICES


###########
# Constants
###########


ENCODING_CHOICES = Choices(
    ('custom', 'CUSTOM', 'Custom'),
    ('dna', 'DNA', 'DNA'),
    ('rna', 'RNA', 'RNA')
)
SEQUENCE_TYPE_CHOICES = Choices(
    ('slr', 'SINGLE_LETTER', 'Single Letter Representation'),
    ('flr', 'FOUR_LETTER', 'Four Letter Representation')
)
MODIFICATION_CHOICES = Choices(
    ('none', 'NONE', 'None'),
    ('standard', 'STANDARD', 'Standard (3/3 nucleotides)'),
    ('ultra', 'ULTRA', 'Ultra (47 nucleotides), for 100mers only'),
)
RNA_BASES = ['A', 'C', 'G', 'U']
DNA_BASES = ['A', 'C', 'G', 'T']
LINKAGES = ['o', 's']
MODIFIERS = ['-', 'm', 'b']
BACKBONES = ['r', 'd', 'm']
IGNORED_CHARACTERS = [' ']  # these get removed from the sequence
MINIMUM_SEQUENCE_LENGTH = 10
ULTRA_MOD_REQUIRED_SEQUENCE_LENGTH = 100

NUCLEOTIDE_ELEMENTS = Choices(
    ('modifier', 'MODIFIER', 'Modifier'),
    ('base', 'BASE', 'Base'),
    ('backbone', 'BACKBONE', 'Backbone'),
    ('linkage', 'LINKAGE', 'Linkage'),
)

SGRNA_SLR_SEQUENCE_SUFFIX = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"
PRODUCT_SLR_SEQUENCE_SUFFIXES = {
    PRODUCT_SLUG_CHOICES.SGRNA_CELL_VALIDATED: SGRNA_SLR_SEQUENCE_SUFFIX,
    PRODUCT_SLUG_CHOICES.SGRNA_EDITED_CELLS: SGRNA_SLR_SEQUENCE_SUFFIX,
    PRODUCT_SLUG_CHOICES.SGRNA_KIT: SGRNA_SLR_SEQUENCE_SUFFIX,
    PRODUCT_SLUG_CHOICES.SGRNA_NOBUFFER: SGRNA_SLR_SEQUENCE_SUFFIX,
    PRODUCT_SLUG_CHOICES.SGRNA_SCREENING_PLATE: SGRNA_SLR_SEQUENCE_SUFFIX,
}


############
# Validators
############


def validate_sequence_type(sequence_type: str) -> str:
    if sequence_type and str(sequence_type).lower() in SEQUENCE_TYPE_CHOICES:
        return str(sequence_type).lower()
    raise SequenceValidationException('sequence_type {} is not supported'.format(sequence_type))


def validate_modification(modification: str, sequence_length: int) -> str:
    if modification:
        if modification.lower() == MODIFICATION_CHOICES.ULTRA and sequence_length != ULTRA_MOD_REQUIRED_SEQUENCE_LENGTH:
            raise SequenceValidationException('Ultra mod is only supported for {}mers, sequence length is {}'.format(
                ULTRA_MOD_REQUIRED_SEQUENCE_LENGTH, sequence_length))
        if modification.lower() in MODIFICATION_CHOICES:
            return modification.lower()
    raise SequenceValidationException('Modification {} is not supported'.format(modification))


def validate_encoding(encoding: str) -> str:
    if encoding:
        if encoding.lower() in ENCODING_CHOICES:
            return str(encoding).lower()
        elif encoding.lower() == 'simple_dna':
            return ENCODING_CHOICES.DNA
        elif encoding.lower() == 'simple_rna':
            return ENCODING_CHOICES.RNA
    raise SequenceValidationException('encoding {} is not supported'.format(encoding))


def infer_encoding(raw_sequence: str) -> str:
    clean_raw_sequence = clean_sequence(raw_sequence)
    if clean_raw_sequence is None or clean_raw_sequence == '':
        return ENCODING_CHOICES.RNA
    rna_base_count = reduce(lambda count, m: count + 1 if m in RNA_BASES else count, clean_raw_sequence.upper(), 0)
    dna_base_count = reduce(lambda count, m: count + 1 if m in DNA_BASES else count, clean_raw_sequence.upper(), 0)
    if dna_base_count == 0 and rna_base_count == 0:
        raise SequenceValidationException('Unable to determine encoding from sequence: {}'.format(clean_raw_sequence))
    elif dna_base_count > rna_base_count:
        return ENCODING_CHOICES.DNA
    elif rna_base_count > dna_base_count:
        return ENCODING_CHOICES.RNA
    else:
        return ENCODING_CHOICES.RNA


def clean_sequence(raw_sequence: str) -> str:
    clean_raw_sequence = str(raw_sequence)
    for c in IGNORED_CHARACTERS:
        clean_raw_sequence = clean_raw_sequence.replace(c, '')
    return clean_raw_sequence


def infer_sequence_type(raw_sequence: str) -> str:
    clean_raw_sequence = clean_sequence(raw_sequence)
    if clean_raw_sequence is None or clean_raw_sequence == '':
        return SEQUENCE_TYPE_CHOICES.SINGLE_LETTER
    modifier_count = reduce(lambda count, m: count + 1 if m in MODIFIERS else count, clean_raw_sequence.lower(), 0)
    if modifier_count > 0:
        return SEQUENCE_TYPE_CHOICES.FOUR_LETTER
    base_count = reduce(lambda count, b: count + 1 if b in set(DNA_BASES + RNA_BASES) else count, clean_raw_sequence.upper(), 0)
    non_base_count = reduce(lambda count, b: count + 1 if b not in set(DNA_BASES + RNA_BASES) else count, clean_raw_sequence.upper(), 0)
    if non_base_count > 0 or base_count == 0:
        raise SequenceValidationException('Unable to determine sequence_type from sequence: {}'.format(clean_raw_sequence))
    return SEQUENCE_TYPE_CHOICES.SINGLE_LETTER


def validate_sequence(raw_sequence: str, encoding: str, sequence_type: str, normalize: bool = True) -> str:
    if raw_sequence is None or raw_sequence.strip() == '':
        return ''
    elif sequence_type == SEQUENCE_TYPE_CHOICES.SINGLE_LETTER:
        return validate_single_letter_sequence(raw_sequence, encoding, normalize)
    elif sequence_type == SEQUENCE_TYPE_CHOICES.FOUR_LETTER:
        return validate_four_letter_sequence(raw_sequence, encoding, normalize)


def validate_single_letter_sequence(raw_sequence: str, encoding: str, normalize: bool = True) -> str:
    validated_bases = {
        ENCODING_CHOICES.CUSTOM: set(RNA_BASES + DNA_BASES),
        ENCODING_CHOICES.DNA: DNA_BASES,
        ENCODING_CHOICES.RNA: RNA_BASES
    }[encoding]
    clean_raw_sequence = clean_sequence(raw_sequence)
    if len(clean_raw_sequence) < MINIMUM_SEQUENCE_LENGTH:
        raise SequenceValidationException('Sequence length is {} and is below minimum of {}'.format(
            len(clean_raw_sequence), MINIMUM_SEQUENCE_LENGTH))
    sequence = ""
    for i, c in enumerate(clean_raw_sequence):
        if c.upper() in validated_bases:
            sequence += c.upper() if normalize else c
        else:
            raise SequenceValidationException("Invalid base '{}' found at position {}".format(c, i))
    return sequence


def validate_four_letter_sequence(raw_sequence: str, encoding: str, normalize: bool = True) -> str:
    validated_bases = {
        ENCODING_CHOICES.CUSTOM: set(RNA_BASES + DNA_BASES),
        ENCODING_CHOICES.DNA: DNA_BASES,
        ENCODING_CHOICES.RNA: RNA_BASES
    }[encoding]
    clean_raw_sequence = clean_sequence(raw_sequence)
    modifier_count = reduce(lambda count, m: count + 1 if m in MODIFIERS else count, clean_raw_sequence.lower(), 0)
    if modifier_count > 0 and clean_raw_sequence[0].lower() not in MODIFIERS:
        raise SequenceValidationException("Modifier is missing at start of sequence")  # common excel error
    sequence = ""
    base_count = 0
    for i, c in enumerate(clean_raw_sequence):
        if i % 4 == 0:
            if c.lower() not in MODIFIERS:
                raise SequenceValidationException("Invalid modifier '{}' found at position {}".format(c, i))
            sequence += c.lower() if normalize else c
        elif i % 4 == 1:
            if c.upper() not in validated_bases:
                raise SequenceValidationException("Invalid base '{}' found at position {}".format(c, i))
            sequence += c.upper() if normalize else c
            base_count += 1
        elif i % 4 == 2:
            if c.lower() not in BACKBONES:
                raise SequenceValidationException("Invalid backbone '{}' found at position {}".format(c, i))
            sequence += c.lower() if normalize else c
        elif i % 4 == 3:
            if c.lower() not in LINKAGES:
                raise SequenceValidationException("Invalid linkage '{}' found at position {}".format(c, i))
            sequence += c.lower() if normalize else c
        if (i + 1) == len(clean_raw_sequence) and i % 4 != 2:
            raise SequenceValidationException("Linkage found at end of sequence]")
    if base_count < MINIMUM_SEQUENCE_LENGTH:
        raise SequenceValidationException('Sequence length is {} and is below minimum of {}'.format(
            base_count, MINIMUM_SEQUENCE_LENGTH))
    return sequence


#######
# Utils
#######


def modify_four_letter_sequence(four_letter_sequence: str, modification_type: str, nucleotide_index_array: list, modification: str) -> str:
    cleaned_sequence = clean_sequence(four_letter_sequence)
    sequence = ''
    nucleotide_index = 1
    for i, c in enumerate(cleaned_sequence):
        if i % 4 == 0:
            if modification_type == NUCLEOTIDE_ELEMENTS.MODIFIER and nucleotide_index in nucleotide_index_array:
                sequence += modification
            else:    
                sequence += c
        elif i % 4 == 1:
            if modification_type == NUCLEOTIDE_ELEMENTS.BASE and nucleotide_index in nucleotide_index_array:
                sequence += modification
            else:    
                sequence += c
        elif i % 4 == 2:
            if modification_type == NUCLEOTIDE_ELEMENTS.BACKBONE and nucleotide_index in nucleotide_index_array:
                sequence += modification
            else:    
                sequence += c
        elif i % 4 == 3:
            if modification_type == NUCLEOTIDE_ELEMENTS.LINKAGE and nucleotide_index in nucleotide_index_array:
                sequence += modification
            else:    
                sequence += c
            nucleotide_index += 1
    return sequence


############
# Main Class
############


class Sequence(object):

    """
    GOALS:
    1) storage of customer sequences in json
    2) loading from json
    3) manipulation and generation of sequence attributes
    4) generation of sequence json representation required by factory systems
    5) validation

    NOTABLE SEQUENCE ATTRIBUTES:
    
    sequence_type: 
    this is the type of the customer/raw sequence
    choices: Single letter representation (SLR), four letter represensation (FLR)
    SLR example: 'GATTACA' (note this is a DNA example, since it includes the letter T)
    FLR example: '-Gdo-Ado-Tdo-Tdo-Ado-Cdo-Ad' (this is the FLR for SLR above)
    Note that the factory requires a sequence in FLR for processing, 
    so all sequences need to be converted to FLR eventually

    encoding:
    the encoding helps determine how to generate the FLR from SLR,
    specifically which backbone to use, ie 'r' for RNA and 'd' for DNA 
    encoding choices are: 'dna', 'rna', 'custom'
    custom denotes sequences that are to be provided as FLR only and do not adhere to a strict format

    raw_sequence: 
    this is the sequence that was provided by the customer.
    RNA sequences provided in SLR are a combination of the following four bases: 'A', 'C', 'G', 'U'
    DNA sequences provided in SLR can also use the base 'T' in addition to the four RNA bases
    sequences provided in FLR format includes bases encoded as follows:
    Modifier + Base + Backbone + Linkage, eg '-Gdo'
    
    modified:
    this determines whether the sequence should be modified/stabilized
    if true, the backbone in the first three and the last three bases is changed to 'm' (from 'r' or 'd')

    product_sequence:
    product can be optionally provided to this class,
    some synthego products append a sequence suffix to the customer sequence, see PRODUCT_SLR_SEQUENCE_SUFFIXES
    in this case, the customer sequence must be provided in SLR, so both sequences can be concatenated
    """

    def __init__(self, 
        raw_sequence: str,
        encoding: str = ENCODING_CHOICES.RNA,
        sequence_type: str = None,
        modification: str = MODIFICATION_CHOICES.NONE,
        label: str = None,
        product: object = None
    ):
        self.sequence_type = validate_sequence_type(sequence_type) if sequence_type else infer_sequence_type(raw_sequence)
        self.encoding = validate_encoding(encoding) if encoding else infer_encoding(raw_sequence)
        self.raw_sequence = validate_sequence(raw_sequence, self.encoding, self.sequence_type)
        self.label = label
        self.product = product
        self.set_product_sequence()
        if self.sequence_type == SEQUENCE_TYPE_CHOICES.SINGLE_LETTER:
            if self.encoding == ENCODING_CHOICES.RNA:
                self.chemistry_object = ChemistrySequence.new(
                    three_letter_representation=ChemistrySequence.rna_to_three_letter(self.raw_sequence + self.product_sequence))
            else:
                self.chemistry_object = ChemistrySequence.new(single_letter_representation=self.raw_sequence + self.product_sequence)
        else:
            self.chemistry_object = ChemistrySequence.new(three_letter_representation=self.raw_sequence or self.product_sequence)
        self.modification = validate_modification(modification, self.sequence_length)

    @property
    def four_letter_representation(self) -> str:
        if self.modified:
            return Sequence.modify(self.chemistry_object.three_letter_representation, self.modification)
        else:
            return self.chemistry_object.three_letter_representation

    @classmethod
    def from_json(cls, data: object, product=None) -> object:
        try:
            return cls(
                raw_sequence=data['customer_sequence'],
                encoding=data['sequence_encoding'],
                sequence_type=data['sequence_type'],
                modified=data['modified'],
                label=data['customer_label'],
                product=product
            )
        except KeyError as e:
            raise SequenceValidationException('Unable to load sequence from json, missing key {} in {}'.format(e, data))

    @property
    def mass(self) -> int:
        return int(self.chemistry_object.mass)

    @property
    def modified(self) -> bool:
        return self.modification != MODIFICATION_CHOICES.NONE

    @classmethod
    def modify(cls, four_letter_sequence: str, modification: str) -> str:
        modified_sequence = validate_four_letter_sequence(four_letter_sequence)
        sequence_length = ChemistrySequence.new(three_letter_representation=modified_sequence).length
        if modification == MODIFICATION_CHOICES.STANDARD:
            """
            Changes the backbone and linkage for the first three and last three bases of an FLR sequence,
            eg, changes
            '-Gro-Aro-Uro-Uro-Aro-Cro-Ar'
            to
            '-Gms-Ams-Ums-Uro-Ams-Cms-Am'
            """
            modified_sequence = modify_four_letter_sequence(
                modified_sequence, 
                NUCLEOTIDE_ELEMENTS.BACKBONE, 
                [1, 2, 3, sequence_length - 2, sequence_length - 1, sequence_length],
                'm'
            )
            modified_sequence = modify_four_letter_sequence(
                modified_sequence, 
                NUCLEOTIDE_ELEMENTS.LINKAGE, 
                [1, 2, 3, sequence_length - 2, sequence_length - 1, sequence_length],
                's'
            )
        elif modification == MODIFICATION_CHOICES.ULTRA:
            """
            Sample result:
            -Ams-Ams-Gms-Uro-Aro-Aro-Aro-Aro-Cro-Cro-Uro-Cro-Uro-Aro-Cro-Aro-Aro-Aro-Uro-Gro-Gro-Uro-Uro-Uro-Uro-Aro-Gro-Aro-Gmo-Cmo-Umo-Amo-Gmo-Amo-Amo-Amo-Umo-Amo-Gmo-Cmo-Aro-Aro-Gro-Uro-Uro-Aro-Aro-Aro-Aro-Uro-Aro-Aro-Gro-Gro-Cro-Uro-Aro-Gro-Uro-Cro-Cro-Gro-Uro-Uro-Aro-Uro-Cro-Aro-Amo-Cmo-Umo-Umo-Gmo-Amo-Amo-Amo-Amo-Amo-Gmo-Umo-Gmo-Gmo-Cmo-Amo-Cmo-Cmo-Gmo-Amo-Gmo-Umo-Cmo-Gmo-Gmo-Umo-Gmo-Cmo-Ums-Ums-Ums-Um            
            """
            if sequence_length != ULTRA_MOD_REQUIRED_SEQUENCE_LENGTH:
                raise SequenceValidationException('Ultra mod is only supported for {}mers, sequence length is {}'.format(
                    ULTRA_MOD_REQUIRED_SEQUENCE_LENGTH, sequence_length))
            modified_sequence = modify_four_letter_sequence(
                modified_sequence, 
                NUCLEOTIDE_ELEMENTS.BACKBONE, 
                [1, 2, 3, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100],
                'm'
            )
            modified_sequence = modify_four_letter_sequence(
                modified_sequence, 
                NUCLEOTIDE_ELEMENTS.LINKAGE, 
                [1, 2, 3, 97, 98, 99, 100],
                's'
            )
        return modified_sequence

    @property
    def sequence_length(self) -> int:
        return self.chemistry_object.length

    @property
    def sequence_length_range(self) -> str:
        return get_sequence_length_range(self.sequence_length)

    def set_product_sequence(self) -> None:
        self.product_sequence = ""
        if self.product is None:
            return
        if self.product.slug in PRODUCT_SLR_SEQUENCE_SUFFIXES:
            if self.sequence_type != SEQUENCE_TYPE_CHOICES.SINGLE_LETTER:
                raise SequenceValidationException('Sequence must be in single letter format for product {}'.format(self.product.slug))
            self.product_sequence = PRODUCT_SLR_SEQUENCE_SUFFIXES[self.product.slug]

    def to_json(self) -> object:
        return {
            'customer_label': self.label,
            'customer_sequence': self.raw_sequence,
            'four_letter_sequence': self.four_letter_representation,
            'modified': self.modified,
            'product_sequence': self.product_sequence,
            'sequence_encoding': self.encoding,
            'sequence_length_range': self.sequence_length_range,
            'sequence_length': self.sequence_length,
            'sequence_type': self.sequence_type,
        }
