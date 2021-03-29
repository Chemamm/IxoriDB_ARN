from django.shortcuts import render, redirect
from .forms import ProteinIDForm, FeaturesForm, ProteinBlastForm, TranscriptIDForm
from transcripts.models import Transcript
from .utils import blastCommandLine, Blast, intersection
from IxoriDB.settings import SUBSITE, UPLOAD_DIR
import os

# Create your views here.


def home_view(request):
    return render(request, "home.html", {})

def getstarted_view(request, *args, **kwargs):
    return render(request, "getstarted.html", {})

def cite_view(request, *args, **kwargs):
    return render(request, "cite.html", {})

def about_view(request, *args, **kwargs):
    return render(request, "about.html", {})

def contact_view(request, *args, **kwargs):
    return render(request, "contact.html", {})

def browse_protein_view(request, *args, **kwargs):
    formID = ProteinIDForm(request.POST or None)
    context = {}
    if request.method == 'POST':
        if formID.is_valid():
            proteins = list(Transcript.objects.filter(hitID=request.POST.get("ProteinID")))
            context["proteins"] = proteins
            return render(request, "browse_proteinID_results.html", context)

    context = {'formID': formID}
    return render(request, "browse_protein.html", context)

def browse_transcript_view(request, *args, **kwargs):
    formID = TranscriptIDForm(request.POST or None)
    context = {}
    if request.method == 'POST':
        if formID.is_valid():
            ID = request.POST.get("TranscriptID")
            return redirect(SUBSITE + '/browse_transcript/%s' % ID)
    context = {'formID': formID}
    return render(request, "browse_transcript.html", context)


def browse_transcript_result_view(request, ID, *args, **kwargs):
    context = {}
    if len(ID.split(".")) > 1:
        ID = ID.split(".")[0]
    transcript = Transcript.objects.get(transcriptID=ID)
    context["transcript"] = transcript
    return render(request, "browse_transcript_results.html", context)


def browse_by_feature_view(request,*args, **kwargs):
    context={}
    formFeatures = FeaturesForm(request.POST or None)
    if request.method == 'POST':
        if formFeatures.is_valid():
            transcripts = list(Transcript.objects.all())
            if request.POST.get("family"):
                lst = list(Transcript.objects.filter(family__icontains=request.POST.get("family")))
                transcripts = intersection(transcripts,lst)
            if request.POST.get("secreted"):
                if request.POST.get("secreted") != "Not apply":
                    lst = list(Transcript.objects.filter(classification__icontains="Putatively Secreted"))
                    transcripts = intersection(transcripts,lst)
            if request.POST.get("transmembrane"):
                if request.POST.get("transmembrane") != "Not apply":
                    lst = list(Transcript.objects.exclude(tm_domain_all__icontains="0"))
                    transcripts = intersection(transcripts,lst)
            if request.POST.get("classification"):
                if request.POST.get("classification") != "Not apply":
                    lst = list(Transcript.objects.filter(classification__icontains=request.POST.get("classification")))
                    transcripts = intersection(transcripts,lst)
            if request.POST.get("eclass"):
                if request.POST.get("eclass") != "Not apply":
                    lst = list(Transcript.objects.filter(eclass__icontains=request.POST.get("eclass")))
                    transcripts = intersection(transcripts,lst)
            if request.POST.get("goterm"):
                lst = list(Transcript.objects.filter(merged__icontains=request.POST.get("goterm")))
                transcripts = intersection(transcripts,lst)
            context["transcripts"] = transcripts
            return render(request,"browse_features_results.html", context)


    context["formFeatures"] = formFeatures
    return render(request, "browse_features.html", context)



def blast_view(request, *args, **kwargs):
    formBlast = ProteinBlastForm(request.POST or None, request.FILES or None)
    context = {}
    if formBlast.is_valid():
        sequence = request.POST.get("sequence")
        if not sequence:
            file_seq = request.FILES["file_seq"]
        else:
            file_seq=False
        typeofseq = request.POST.get("typeofsequence")
        evalue = request.POST.get("evalue")
        word_size = request.POST.get("word_size")
        max_target_seqs = request.POST.get("max_target_seqs")

        if file_seq:
            sequence=file_seq.read().decode('utf-8')

        work, fasta_check = blastCommandLine(sequence, typeofseq, evalue, word_size, max_target_seqs)
        context["work"] = work
        if fasta_check:
            blast = Blast(UPLOAD_DIR + "tmp/%i.blast" %work)
            print(blast.hits)
            context["hits"] = blast.hits
            return render(request, "blast_results.html", context)
        else:
            return render(request, "blast_error.html", context)
    context = {'formBlast':formBlast}
    print(formBlast.is_valid())
    return render(request, "blast.html", context)

def download_view(request, *args, **kwargs):
    return render(request, "download.html", {})
