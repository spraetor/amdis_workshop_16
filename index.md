---
layout: page
title: AMDiS Workshop 2016
hidden: 1
---

## Einführung in die Finite-Elemente Software AMDiS

Der AMDiS Workshop bietet eine Einführung in die Finite-Elemente Software AMDiS, 
die am Institut für Wissenschaftliches Rechnen entwickelt wurde und aktiv 
weiterentwickelt wird. Ziel ist es, dass die Teilnehmer den ersten Umgang mit der 
Software erlernen und dann in der Lage sind, einfache aber auch komplexere 
partielle Differentialgleichungen zu lösen. 

Stationäre oder instationäre Probleme, in 1d, 2d, 3d oder auf gekrümmten Oberflächen, 
mit adaptiven Gittern, sequentiell oder parallel - die Teilnehmer sollen einen 
Überblick über verschiedene Techniken erhalten und einen Ausblick auf fortgeschrittene 
Themen bekommen. 

### Modul

Der Workshop ist als **Math Ma WIA Modul** für Studenten anrechnungsfähig. Einen 
Seminarschein kann man erhalten durch Teilnahme und Bearbeitung eines 
Programmierprojekts, mit Vorstellung der Ergebnisse am Ende der Woche. Die Projekte 
werden in der Woche durchgeführt und betreut und bauen auf den Grundlagenvorlesungen 
des Workshops auf.

### Zeit und Ort

08.08.2016 - 12.08.2016, WIL/B221

### Workshop Material

Das Material is Verfügbar auf GitHub: 
[spraetor/amdis_workshop_16](https://github.com/spraetor/amdis_workshop_16) und 
kann mittels `git clone REPOSITORY` heruntergeladen werden.

### Anmeldung

Wir bitten um eine kurze Anmeldung an 
[simon.praetorius@tu-dresden.de](mailto:simon.praetorius@tu-dresden.de) oder über

<form action="//formspree.io/s.praetorius@googlemail.com" method="POST">
<input type="text" name="name" placeholder="Name" /><br />
<input type="email" name="_replyto" placeholder="E-Mail" /><br />
<button type="submit">Anmelden</button><br />
<input type="hidden" name="_next" value="{{ site.url }}/anmeldung" />
<input type="hidden" name="_subject" value="AMDiS-Workshop Anmeldung" />
<input type="text" name="_gotcha" style="display:none" />
</form> 

### Literaturhinweise

- ALBERT 1.0 Documentation: [Download](http://www.alberta-fem.de/download.html) 
  (Alfred Schmidt and Kunibert G. Siebert)
- C++ kurz & gut, Kyle Loudon und Rainer Grimm, 2014

### Software

- [AMDiS](https://fusionforge.zih.tu-dresden.de/projects/amdis):
  Adaptive Multi-Dimensional Simulations ([Download](http://www.math.tu-dresden.de/~spraetor/AMDIS-0.9.3144-Linux.deb))
- [MTL4](http://www.simunova.com/de/node/65): 
  Matrix Template Library, elegant and efficient linear algebra expressions in C++
- [PETSc](http://www.mcs.anl.gov/petsc): 
  Portable, Extensible Toolkit for Scientific Computation
