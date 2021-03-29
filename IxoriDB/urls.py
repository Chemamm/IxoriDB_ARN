"""IxoriDB URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path
from django.conf.urls import include,url
from django.conf.urls.static import static
from django.conf import settings
from pages.views import home_view, browse_protein_view, blast_view, browse_transcript_view, browse_transcript_result_view, \
    browse_by_feature_view, getstarted_view, about_view, contact_view, download_view, cite_view
from IxoriDB.settings import SUBSITE
from pages import views

urlpatterns = [
    path(SUBSITE + 'admin/', admin.site.urls),
    url(r'^$', views.home_view, name='home'),
    url(r'^get_started$', views.getstarted_view, name='getstarted'),
    url(r'^about$', views.about_view, name='about'),
    url(r'^contact$', views.contact_view, name='contact'),
    url(r'^download$', views.download_view, name='download'),
    url(r'^cite$', views.cite_view, name='cite'),
#    path(SUBSITE + 'browse_transcript/<str:ID>', browse_transcript_result_view, name="btranscriptresult"),
    url(r'^browse_transcript/(?P<ID>[\w-]+)', views.browse_transcript_result_view, name='btranscriptresult'),
    url(r'^browse_transcript/$', views.browse_transcript_view, name='btranscript'),
    url(r'^browse_protein$', views.browse_protein_view, name='bprotein'),
    url(r'^browse_by_features$', views.browse_by_feature_view, name='bfeatures'),
    url(r'^blast$', views.blast_view, name='blast'),


]

#urlpatterns = [
#    path(SUBSITE + 'admin/', admin.site.urls),
#    path(SUBSITE + '', home_view, name="home"),
#    path(SUBSITE + 'get_started/', getstarted_view, name="getstarted"),
#    path(SUBSITE + 'about/', about_view, name="about"),
#    path(SUBSITE + 'contact/', contact_view, name="contact"),
#    path(SUBSITE + 'download/', download_view, name="download"),
#    path(SUBSITE + 'cite/', cite_view, name="cite"),
#    path(SUBSITE + 'browse_protein/', browse_protein_view, name="bprotein"),
#    path(SUBSITE + 'browse_transcript/', browse_transcript_view, name="btranscript"),
#    path(SUBSITE + 'browse_transcript/<str:ID>', browse_transcript_result_view, name="btranscriptresult"),
#    path(SUBSITE + 'browse_by_features/', browse_by_feature_view, name="bfeatures"),
#    path(SUBSITE + 'blast/',blast_view, name="blast"),
#]

