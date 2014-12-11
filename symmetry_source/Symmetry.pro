TEMPLATE = app
CONFIG += console
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++0x

INCLUDEPATH += /usr/include/openbabel-2.0

LIBS += -lopenbabel

SOURCES += main.cpp \
    mapperfunctor.cpp \
    symmetry.cpp

HEADERS += \
    mapperfunctor.h \
    symmetry.h

OTHER_FILES += \
    Readme.txt

