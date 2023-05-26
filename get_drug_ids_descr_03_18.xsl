<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">



<xsl:output
method="text" omit-xml-declaration="yes" indent="yes"/> 



<xsl:template match="/">
<xsl:apply-templates select="drugbank"/>
</xsl:template>


<xsl:template match="drugbank">
<xsl:apply-templates select="drug[@type='small molecule']"/>
</xsl:template>


<xsl:template match="drug[@type='small molecule']">
<xsl:value-of select="concat(drugbank-id[@primary='true'],'|',name)"/>
<xsl:value-of select="concat('|',groups/group)"/>
<xsl:for-each select="calculated-properties/property">
<xsl:if test="kind[text()='MDDR-Like Rule']">
<xsl:value-of select="concat('|',value)"/>
</xsl:if>
<xsl:if test="kind[text()='Rule of Five']">
<xsl:value-of select="concat('|',value)"/>
</xsl:if>
<xsl:if test="kind[text()='SMILES']">
<xsl:value-of select="concat('|',value)"/>
</xsl:if>
</xsl:for-each>
<xsl:text>
</xsl:text>
</xsl:template>


</xsl:stylesheet>
