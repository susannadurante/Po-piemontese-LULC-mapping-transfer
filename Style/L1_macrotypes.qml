<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34" styleCategories="Symbology">
  <pipe>
    <rasterrenderer type="paletted" band="1" opacity="1" alphaBand="-1" nodataColor="">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
      </minMaxOrigin>
      <colorPalette>
        <paletteEntry value="1" color="#1565c0" alpha="255" label="1 - Acque e zone umide"/>
        <paletteEntry value="2" color="#c9b896" alpha="255" label="2 - Substrato fluviale"/>
        <paletteEntry value="3" color="#b8d08a" alpha="255" label="3 - Vegetazione erbacea e arbustiva"/>
        <paletteEntry value="4" color="#2e7d32" alpha="255" label="4 - Copertura arborea"/>
        <paletteEntry value="5" color="#f9a825" alpha="255" label="5 - Superfici coltivate"/>
        <paletteEntry value="6" color="#c62828" alpha="255" label="6 - Superfici impermeabilizzate"/>
        <paletteEntry value="7" color="#616161" alpha="255" label="7 - Aree estrattive"/>
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation saturation="0" grayscaleMode="0" colorizeStrength="100" colorizeOn="0"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
