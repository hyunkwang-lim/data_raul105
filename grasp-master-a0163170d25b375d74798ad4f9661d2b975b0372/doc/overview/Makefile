################################################################
#                                                              #
#         Makefile for building the GRASP documentation        #
#                                                              #
#      Fabrice Ducos, Laboratoire d'Optique Atmospherique      #
#                                                              #
#                fabrice.ducos@univ-lille1.fr                  #
#                                                              #
#                     18 November 2013                         #
#                                                              #
# Requirements:                                                #
# For using this Makefile in order to build the documentation, #
# one needs:                                                   #
# - an XSLT processor (xsltproc by default, or alternate tools #
# Saxon or Xalan)                                              #
# - a Formatting Objects processor (currently Apache fop).     #
# - The Apache Xerces XML parser (for use with Saxon or Xalan) #
#                                                              #
# These tools are not provided with the project and must be    #
# installed on the system (in the PATH for xsltproc and fop,   #
# or in the CLASSPATH for the java-based tools Saxon and       #
# Xalan).                                                      #
#                                                              #
#                                                              #
################################################################

# Saxon and Xalan will look for the project's CatalogManager.properties file
# in the current directory
export CLASSPATH:=$(CLASSPATH):.

XERCES_FLAGS_1= \
	-Djavax.xml.parsers.DocumentBuilderFactory=org.apache.xerces.jaxp.DocumentBuilderFactoryImpl \
	-Djavax.xml.parsers.SAXParserFactory=org.apache.xerces.jaxp.SAXParserFactoryImpl
XERCES_FLAGS_2= \
	-Dorg.apache.xerces.xni.parser.XMLParserConfiguration=org.apache.xerces.parsers.XIncludeParserConfiguration

RESOLVER_NAMESPACE=org.apache.xml.resolver.tools
SAXON_CATALOG_FLAGS= \
	-x $(RESOLVER_NAMESPACE).ResolvingXMLReader \
	-y $(RESOLVER_NAMESPACE).ResolvingXMLReader \
	-r $(RESOLVER_NAMESPACE).CatalogResolver \
	-u

XALAN_CATALOG_FLAGS= \
	-ENTITYRESOLVER $(RESOLVER_NAMESPACE).CatalogResolver \
	-URIRESOLVER $(RESOLVER_NAMESPACE).CatalogResolver

# for xsltproc
export XML_CATALOG_FILES=catalog.xml
export XML_DEBUG_CATALOG=1

#lang=en
#srcdir=src/$(lang) # for multilanguaje system. For the moment, just english is available
srcdir=src
docdir=doc

entities=$(srcdir)/common_entities.dtd
main=$(srcdir)/main.dbk
chapters=$(srcdir)/preface.dbk \
         $(srcdir)/chap01.dbk \
         $(srcdir)/chap02.dbk \
         $(srcdir)/chap03.dbk \
         $(srcdir)/chap04.dbk \
	 $(srcdir)/biblio.dbk \
	 $(srcdir)/glossary.dbk

htdocs=$(PWD)
figures=figures
website=$(htdocs)/website

XSLTPROC_COMMON_FLAGS=--xinclude --stringparam section.autolabel 1 --stringparam  section.label.includes.component.label 1
XSLTPROC=xsltproc

SAXON_COMMON_FLAGS=section.autolabel=1 verbosity=4
SAXON=java $(XERCES_FLAGS_1) $(XERCES_FLAGS_2) com.icl.saxon.StyleSheet $(SAXON_CATALOG_FLAGS)

XALAN_COMMON_FLAGS=-param section.autolabel 1
XALAN=java $(XERCES_FLAGS_2) org.apache.xalan.xslt.Process $(XALAN_CATALOG_FLAGS)

.PHONY: all web pdf single-html chunk-html clean

all: pdf single-html chunk-html

pdf: $(docdir)/docgrasp.pdf

grasp-open: chunk-grasp-open

$(docdir)/docgrasp.pdf: $(docdir)/docgrasp.fo
	fop $< $@

$(docdir)/docgrasp.fo: $(main) $(docdir) $(entities) $(chapters)
	$(XSLTPROC) --output $@ $(XSLTPROC_COMMON_FLAGS) --stringparam fop1.extensions 1 --stringparam paper.type A4 fo $<
	@#$(SAXON) -o $@ $< fo $(SAXON_COMMON_FLAGS) fop1.extensions=1 paper.type=A4
	@#$(XALAN) -out $@ -in $< -xsl fo $(XALAN_COMMON_FLAGS) -param fop1.extensions 1 -param paper.type A4

single-html: $(main) $(docdir) $(entities) $(chapters)
	$(XSLTPROC)  $(XSLTPROC_COMMON_FLAGS) --output $(docdir)/index.html single-html $<
	@#$(SAXON) -o $(docdir)/index.html $< single-html $(SAXON_COMMON_FLAGS)
	@#$(XALAN) -out $(docdir)/index.html -in $< -xsl single-html $(XALAN_COMMON_FLAGS)

chunk-html: $(main) $(entities) $(chapters) $(website)
	$(XSLTPROC) $(XSLTPROC_COMMON_FLAGS) --stringparam base.dir $(website)/ chunk-html $<
	@#$(SAXON) $< chunk-html $(SAXON_COMMON_FLAGS) base.dir=$(website)/
	@#$(XALAN) -in $< -xsl chunk-html $(XALAN_COMMON_FLAGS) -param base.dir $(website)/

chunk-grasp-open: $(main) $(entities) $(chapters) $(website)
	$(XSLTPROC) $(XSLTPROC_COMMON_FLAGS) --stringparam chunker.output.encoding UTF-8 --stringparam html.ext ".php" --stringparam base.dir $(website)/ grasp-open $<
	@#$(SAXON) $< chunk-html $(SAXON_COMMON_FLAGS) base.dir=$(website)/
	@#$(XALAN) -in $< -xsl chunk-html $(XALAN_COMMON_FLAGS) -param base.dir $(website)/
	cp $(docdir)/docgrasp.pdf $(website)/docgrasp.pdf
	cp $(docdir)/../docbook-xsl-ns/grasp-open/docbook-custom-style.css $(website)/docbook-custom-style.css
#	rm $(website)/figures $(website)/images
#	cp -r figures $(website)/figures
#	cp -r docbook-xsl-ns/images $(website)/images	

$(docdir):
	mkdir $(docdir)
	cd $(docdir) && ln -s ../$(figures) . && ln -s ../docbook-xsl-ns/images .

$(website): $(docdir)/docgrasp.pdf
	-mkdir $(website)
	-cd $(website) && ln -s ../$(figures) . && ln -s ../docbook-xsl-ns/images .
	cp $(docdir)/docgrasp.pdf $(website)

clean:
	-rm -f *~ src/*~
	-rm -rf $(docdir)
	-rm -rf $(website)

