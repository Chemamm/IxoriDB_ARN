from django.db import models


# Create your models here.


class Transcript(models.Model):
    Sequence = models.TextField()
    CDS_ID = models.CharField(max_length=100)
    ArachnidaTopHit = models.CharField(blank=True, max_length=100, null=True)
    ArachnidaIdentity = models.DecimalField(blank=True, max_digits=10, decimal_places=3, null=True)
    ArachnidaCoverage = models.IntegerField(blank=True, null=True)
    ArachnidaBitScore = models.DecimalField(blank=True, max_digits=10, decimal_places=1, null=True)
    ArachnidaDescription = models.TextField(blank=True, null=True)
    SprotTopHit = models.CharField(blank=True, max_length=100, null=True)
    SprotIdentity = models.DecimalField(blank=True, max_digits=10, decimal_places=3, null=True)
    SprotCoverage = models.IntegerField(blank=True, null=True)
    SprotBitScore = models.DecimalField(blank=True, max_digits=10, decimal_places=1, null=True)
    SprotDescription = models.TextField(blank=True, null=True)
    Uniref90TopHit = models.CharField(blank=True, max_length=100, null=True)
    Uniref90Identity = models.DecimalField(blank=True, max_digits=10, decimal_places=3, null=True)
    Uniref90Coverage = models.IntegerField(blank=True, null=True)
    Uniref90BitScore = models.DecimalField(blank=True, max_digits=10, decimal_places=1, null=True)
    Uniref90Description = models.TextField(blank=True, null=True)
    TSFTopHit = models.CharField(blank=True, max_length=100, null=True)
    TSFIdentity = models.DecimalField(blank=True, max_digits=10, decimal_places=3, null=True)
    TSFCoverage = models.IntegerField(blank=True, null=True)
    TSFBitScore = models.DecimalField(blank=True, max_digits=10, decimal_places=1, null=True)
    TSFDescription = models.TextField(blank=True, null=True)
    BestHit = models.CharField(max_length=100, null=True)
    BestIdentity = models.DecimalField(max_digits=10, decimal_places=3)
    BestCoverage = models.IntegerField()
    BestBitScore = models.DecimalField(max_digits=10, decimal_places=1)
    BestDescription = models.TextField()
    Secreted = models.BooleanField()
    Transmembrane = models.BooleanField()
    


class ProteinTable:
    class Meta:
        model = Transcript

