Zusammenfassung
+++++++++++++++

Die Schnittstelle zwischen Kraftwerk und geologischem Speicher stellt eine flexible Methode zum Koppeln einer Speichersimulation mit dem Betrieb von Kraftwerksanlagen. Der Speicherbetrieb wird mit Hilfe von Entscheidungsregeln gesteuert und übergibt die Schnittstellenparameter zwischen den Modellen automatisch. Der modulare Aufbau der Softwareschnittstelle erlaubt es, zwischen unterschiedlichen Simulationsmodellen für Speicher und Kraftwerk auszuwählen.

Abstract
++++++++

The Power Plant Geostorage Interface provides a flexible method to couple the simulation of a geological storage to the operation of a power plant. The interface provides a operation control logic and exchanges the interface parameters between the sotrage and the power plant model automatically. Due to a modular architecture, it is possible to choose between different simulation models for the power plant and the geological storage.

ROADMAP
+++++++

Kraftwerksanlagen
-----------------

* Einzelne Kraftwerksanlagen ansteuern, die parallel geschaltet sind.
* Den Massenstrom auf Grundlage der Eintrittstemperatur unter Vorgabe der Austrittstemperatur und des Wärmestroms berechnen.
* Die Massenströme der Einzelanlagen aufaddieren und (inkl. Temperatur) zur Schnittstelle schicken.

Geospeicher
-----------

* Speichersimulation unter Vorgabe der Speichereintrittstemperatur und des Massenstroms (oder Wärmestroms).

Wärmetauscher
-------------
    
* Vorgabe der Wärmeübertragerfläche aus einer Anlagenauslegung (Auslegung auf Wärmestrom X MW bei einer Grädigkeit Y K)
* Berechnung von Speichereintrittstemperatur (und Speichermassenstrom) bei Vorgabe der Netzseitigen Temperaturen sowie des Wärmestroms am Wärmetauscher und unter Annahme der Austrittstemperatur aus dem Wärmespeicher.

Kopplung
--------

Prinzip
^^^^^^^

1 Einlesen/Initialisieren der Modelle, Start der zeitschrittweisen Simulation
2 Feststellen, ob Wärme ein- oder ausgespeichert (oder weder ein- noch ausgespeichert) wird
3 Vorgabe der Wärmeströme für die einzelnen Kraftwerksanlagen und Berechnung des gesamten Massenstroms aus den Anlagen in das Netz sowie in den Speicher
4 Berechnung der Schnittstellenparameter am Wärmeübertrager
5 Speichersimulation zur Bestimmung der tatsächlichen Temperatur im Speicherrücklauf (die Temperatur des Mediums, das aus dem Speicher strömt)
6 Prüfen, ob Annahme am Wärmeübertrager mit Berechnung aus 5 übereinstimmt.
6.1 Stimmen überein -> 7
6.2 Stimmen nicht überein, anpassen der Annahme für 3/4 und im selben Zeitschritt zu 3.
7 Schreiben der Schnittstellenparameter in outputfile und nächster Zeitschritt (2).

Steuerungseingriffe
^^^^^^^^^^^^^^^^^^^

* Limitationen in der Speicherrate.
* Nicht-Erreichbarkeit der gewünschten Netzvorlauftemperatur, da Temperatur im Speicher zu gering (=Speicher leer).
* Keine Einspeicherung möglich, da Speicheraustrittstemperatur zu hoch (=Speicher voll).
