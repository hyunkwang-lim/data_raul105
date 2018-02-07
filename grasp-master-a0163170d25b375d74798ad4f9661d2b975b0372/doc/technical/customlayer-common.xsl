<?xml version='1.0'?> 
<xsl:stylesheet  
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
    xmlns:fo="http://www.w3.org/1999/XSL/Format"
    version="1.0"> 

<!-- my customizations (Fabrice Ducos, 18 november 2013) -->

<!-- ref. http://ds9a.nl/docbook/fancy-docbook.html -->
<!-- Turn on admonition graphics. -->
  <xsl:param name="admon.graphics" select="'1'"/>
  <xsl:param name="admon.graphics.path">images/</xsl:param>

<xsl:param name="callout.graphics" select="'1'"></xsl:param>
<xsl:param name="callout.graphics.path">images/callouts/</xsl:param>

<!-- ref. http://www.sagehill.net/docbookxsl/FittingText.html -->
<xsl:attribute-set name="monospace.verbatim.properties">
  <xsl:attribute name="font-size">
    <xsl:choose>
      <xsl:when test="processing-instruction('db-font-size')"><xsl:value-of
           select="processing-instruction('db-font-size')"/></xsl:when>
      <xsl:otherwise>inherit</xsl:otherwise>
    </xsl:choose>
  </xsl:attribute>
</xsl:attribute-set>

</xsl:stylesheet>
