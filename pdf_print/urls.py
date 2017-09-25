from django.conf.urls import url

from . import views

urlpatterns = [
    url(r'^$', views.cover, name='cover'),
    url(r'^static_out/', views.static_output, name='static_output'),
    url(r'^static_in/', views.static_input, name='static_input')
]

