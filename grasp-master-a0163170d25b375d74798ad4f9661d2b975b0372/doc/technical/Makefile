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

lang=en
mddir=markdown
docdir=doc
docbookdir=docbook

main=$(docbookdir)/index.xml

XSLTPROC_COMMON_FLAGS=--xinclude --stringparam section.autolabel 1
XSLTPROC=xsltproc

SAXON_COMMON_FLAGS=section.autolabel=1 verbosity=4
SAXON=java $(XERCES_FLAGS_1) $(XERCES_FLAGS_2) com.icl.saxon.StyleSheet $(SAXON_CATALOG_FLAGS)

XALAN_COMMON_FLAGS=-param section.autolabel 1
XALAN=java $(XERCES_FLAGS_2) org.apache.xalan.xslt.Process $(XALAN_CATALOG_FLAGS)

htdocs=$(PWD)
figures=figures
website=$(htdocs)/website
webhost='loawww.univ-lille1.fr'
public_dir="/Users/ducos/Sites/public/grasp-tech-doc"

all: prepare-doxygen.sh Doxyfile_template markdown
	./prepare-doxygen.sh html
	doxygen 

doc: prepare-doxygen.sh Doxyfile_template markdown
	./prepare-doxygen.sh php
	doxygen 

public: $(website)
	ssh $(webhost) rm -rf $(public_dir) && scp -r $(website) $(webhost):$(public_dir)

pdf: doc $(docdir)/docgrasp.pdf

$(docdir)/docgrasp.pdf: $(docdir)/docgrasp.fo
	fop $< $@

$(docdir)/docgrasp.fo: $(main)
	$(XSLTPROC) --output $@ $(XSLTPROC_COMMON_FLAGS) --stringparam fop1.extensions 1 --stringparam paper.type A4 fo $<
	@#$(SAXON) -o $@ $< fo $(SAXON_COMMON_FLAGS) fop1.extensions=1 paper.type=A4
	@#$(XALAN) -out $@ -in $< -xsl fo $(XALAN_COMMON_FLAGS) -param fop1.extensions 1 -param paper.type A4

clean-docbook:
	./cleanup.sh

single-html: $(main)
	$(XSLTPROC)  $(XSLTPROC_COMMON_FLAGS) --output $(docdir)/index.html single-html $<
	@#$(SAXON) -o $(docdir)/index.html $< single-html $(SAXON_COMMON_FLAGS)
	@#$(XALAN) -out $(docdir)/index.html -in $< -xsl single-html $(XALAN_COMMON_FLAGS)

grasp-open: doc

clean:
	-rm -rf website latex html
	-rm -f core core.*
	-rm -f *~
	-rm -f *.tmp
	-rm -rf markdown-doxygen 
	-rm -rf docbook
	-rm -rf doc
	-rm -f temp.xml
	-rm Doxyfile Doxyfile-e
