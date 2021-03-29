# Generated by Django 3.1.4 on 2021-01-27 12:19

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Transcript',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('Sequence', models.TextField()),
                ('CDS_ID', models.CharField(max_length=100)),
                ('ArachnidaTopHit', models.CharField(blank=True, max_length=100)),
                ('ArachnidaIdentity', models.DecimalField(blank=True, decimal_places=3, max_digits=10)),
                ('ArachnidaCoverage', models.IntegerField(blank=True)),
                ('ArachnidaBitScore', models.DecimalField(blank=True, decimal_places=1, max_digits=10)),
                ('ArachnidaDescription', models.TextField(blank=True)),
                ('SprotTopHit', models.CharField(blank=True, max_length=100)),
                ('SprotIdentity', models.DecimalField(blank=True, decimal_places=3, max_digits=10)),
                ('SprotCoverage', models.IntegerField(blank=True)),
                ('SprotBitScore', models.DecimalField(blank=True, decimal_places=1, max_digits=10)),
                ('SprotDescription', models.TextField(blank=True)),
                ('Uniref90TopHit', models.CharField(blank=True, max_length=100)),
                ('Uniref90Identity', models.DecimalField(blank=True, decimal_places=3, max_digits=10)),
                ('Uniref90Coverage', models.IntegerField(blank=True)),
                ('Uniref90BitScore', models.DecimalField(blank=True, decimal_places=1, max_digits=10)),
                ('Uniref90Description', models.TextField(blank=True)),
                ('TSFTopHit', models.CharField(blank=True, max_length=100)),
                ('TSFIdentity', models.DecimalField(blank=True, decimal_places=3, max_digits=10)),
                ('TSFCoverage', models.IntegerField(blank=True)),
                ('TSFBitScore', models.DecimalField(blank=True, decimal_places=1, max_digits=10)),
                ('TSFDescription', models.TextField(blank=True)),
                ('BestHit', models.CharField(max_length=100)),
                ('BestIdentity', models.DecimalField(decimal_places=3, max_digits=10)),
                ('BestCoverage', models.IntegerField()),
                ('BestBitScore', models.DecimalField(decimal_places=1, max_digits=10)),
                ('BestDescription', models.TextField()),
                ('Secreted', models.BooleanField()),
                ('Transmembrane', models.BooleanField()),
            ],
        ),
    ]
