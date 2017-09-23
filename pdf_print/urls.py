from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.cover, name='cover'),
    url(r'^pic_output/', views.pic_output)
]

