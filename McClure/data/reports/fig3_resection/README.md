# 2010 avalanche outline from Figure 3 (photo resection)

`McClure/data/release_area/LoneTree_2010_fig3_resection.geojson` holds the 22 Feb 2010 avalanche extent and its release part,
taken from Figure 3 of the CAIC/CDOT incident report (`../201002022_LoneTree Bowl Incident.pdf`).

Method (2026-11-05):
1. Camera pose and focal length fitted to the photo by matching CO 133 (OpenStreetMap), draped on the USGS 1/3 arc-second
   DEM, to the roads traced in the photo (`cam.py`; traced pixels are in the script). The road-only fit put the start zone
   off: the gully in a render of the DEM sat about 60 px from the photo's gully. The final fit adds the centre line of the
   outline's track matched to the gully bottom (from the model flow) between the start zone and the road.
2. The black outline was traced by hand (`extract.py`, `POLY`), every pixel inside it cast onto the DEM, and the hits
   polygonised at 5 m. "Release" = the part above the neck of the outline in the photo (pixel row < 347).
3. Checks: CO 133 and the Crystal River project onto the photo (`check_overlay.png`); the outline spans 2,865-2,353 m
   against the report's 9,400 ft (2,865 m) crown and 7,700 ft (2,347 m) runout; release aspect 53 deg (report: north aspect).

Accuracy: roads 6 px (about 10 m), track 14 px. The start zone is 1 km beyond the roads, so its position is good to
a few tens of metres at best, and the outline in the report was itself hand drawn. The release part (about 2.3 ha) is
larger than "20% of the start zone"; it may include the upper track.

Files: `cam.py` (camera model, DEM lookup, road fit) and `extract.py` (outline ray casting, polygons, stats) are helper
modules; the fitted camera is in the geojson properties. `usgs10.tif` is the USGS 1/3 arc-second DEM clip and
`co133_utm.json` the OSM road lines in EPSG:26913.

Later check (2026-11-05): Ron drew the 2010 release by overlaying Figure 3 on Google Earth
(`data/release_area/LoneTree_2010_release.kml`, 0.53 ha, 23% of the start zone). It lies about 100 m south of this
resection's release part, and model runs from it follow the resection's track only below the gully junction. The upper
part of the resection outline (release and upper track) is therefore off; its lower track and runout agree.
