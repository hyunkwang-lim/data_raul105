<?xml version='1.0'?> 
<xsl:stylesheet  
    xmlns:dbk="http://docbook.org/ns/docbook"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
    xmlns:fo="http://www.w3.org/1999/XSL/Format"
    version="1.0"> 

<xsl:import href="docbook-xsl-ns/fo/docbook.xsl"/> 
<xsl:import href="customlayer-common.xsl"/>
<xsl:param name="admon.graphics.extension">.svg</xsl:param>
<xsl:param name="callout.graphics.extension">.svg</xsl:param>

<!-- My Table Mods. Ref: https://lists.oasis-open.org/archives/docbook-apps/201207/msg00007.html -->
<!-- http://www.oxygenxml.com/forum/topic7848.html -->

<!-- what is below is not validated yet. Do not use. -->
<!--

<xsl:template match="dbk:table[@role='customtable']|dbk:informaltable[@role='customtable']" mode="htmlTable">
  <fo:table background-color="yellow" font-size="10" color="blue">
</xsl:template>
-->

<!-- This works but I can't figure out how to use it locally -->
<xsl:param name="met.table.font.size">0.6</xsl:param>
<xsl:param name="met.table.head.font.size">0.6</xsl:param>

<!-- Set table body font size and alignment -->
<xsl:attribute-set name="table.properties">
  <xsl:attribute
      name="keep-together.within-column">auto</xsl:attribute>
  <xsl:attribute name="font-size">
    <xsl:value-of select="$body.font.master *
			  $met.table.font.size"/>
    <xsl:text>pt</xsl:text>
  </xsl:attribute>
</xsl:attribute-set>

<!-- Set table header font size -->
<xsl:template name="table.row.properties">
  <xsl:if test="ancestor::thead">
    <xsl:attribute name="font-weight">bold</xsl:attribute>
    <xsl:attribute name="color">#FFFFFF</xsl:attribute>
    <xsl:attribute name="background-color">#000000</xsl:attribute>
    <xsl:attribute name="font-size">
      <xsl:value-of
	  select="$body.font.master * $met.table.head.font.size" />
      <xsl:text>pt</xsl:text>
    </xsl:attribute>
  </xsl:if>
</xsl:template>

</xsl:stylesheet>
