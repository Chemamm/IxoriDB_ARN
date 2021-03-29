from django import forms
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, Row, Column, ButtonHolder, HTML
from crispy_forms.bootstrap import Accordion, AccordionGroup

# Create your forms here.


class ProteinIDForm(forms.Form):
    ProteinID = forms.CharField(label="Protein ID:")


class TranscriptIDForm(forms.Form):
    TranscriptID = forms.CharField(label="Transcript ID:")


class FeaturesForm(forms.Form):
    eCHOICES = (
        ("Not apply", "Not apply"),
        ("MG specific", "MG specific"),
        ("SG specific", "SG specific"),
        ("Neutral", "Neutral"),
        )
    CHOICES = (
        ("Not apply", "Not apply"),
        ("Putatively Secreted, Annotated", "Putatively Secreted, Annotated"),
        ("Putatively Non-Secreted, Annotated", "Putatively Non-Secreted, Annotated"),
        ("Putatively Secreted, not Annotated", "Putatively Secreted, not Annotated"),
        ("Putatively Non-Secreted, not Annotated", "Putatively Non-Secreted, not Annotated"),
    )
    sCHOICES = (
        ("Not apply", "Not apply"),
        ("Secreted", "Secreted"),
        ("Non-Secreted", "Non-Secreted"),
    )
    tCHOICES = (
        ("Not apply", "Not apply"),
        ("Transmembrane domains", "Transmembrane domains"),
        ("Not transmembrane domains", "Not transmembrane domains"),
    )
    family = forms.CharField(label="Browse by family:", required=False)
    secreted = forms.ChoiceField(label="Secreted", choices=sCHOICES, required=False)
    transmembrane = forms.ChoiceField(label="Transmembrane Domain", choices=tCHOICES, required=False)
    eclass = forms.ChoiceField(choices=eCHOICES, label="Tissue specificity:", required=False)
    classification = forms.ChoiceField(choices=CHOICES, label="Class:", required=False)
    goterm = forms.CharField(label="Browse by GO term:", required=False)


class ProteinBlastForm(forms.Form):
    sequence = forms.CharField(label="Introduce sequence/s in Fasta format:", widget=forms.Textarea, required=False)
    file_seq = forms.FileField(label="Or select a Fasta file:", required=False)
    CHOICES = (("Peptide", "Peptide"), ("Nucleotide", "Nucleotide"))
    typeofsequence = forms.ChoiceField(choices=CHOICES, label="Select the type of sequence:")
    evalue = forms.DecimalField(label="E-value cut-off", initial=1e-4)
    word_size = forms.IntegerField(label="Word size", initial=3)
    max_target_seqs = forms.IntegerField(label="Max targets per seq", initial=500)


    def clean(self):
        data = self.cleaned_data

        if data.get('sequence', None) and data.get('file_seq', None):
            raise forms.ValidationError('Provide either a Fasta formatted sequence or a Fasta file')
        elif data.get('sequence', None) or data.get('file_seq', None):
            return data
        else:
            raise forms.ValidationError('Provide either a Fasta formatted sequence or a Fasta file')
