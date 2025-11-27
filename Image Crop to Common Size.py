from pathlib import Path
from PIL import Image, ImageFile

# swith of limits for large images (only when certain about folder origin)
Image.MAX_IMAGE_PIXELS = None
ImageFile.LOAD_TRUNCATED_IMAGES = True

in_dir  = Path("./2006 - Luftbilder - Tiff")
out_dir = Path("./2006 - Luftbilder - Cropped")
out_dir.mkdir(exist_ok=True)

# obtain image size
sizes = []
for img_path in in_dir.glob("*.tif"):
    try:
        with Image.open(img_path) as img:
            sizes.append((img_path, img.width, img.height))
    except Exception as e:
        print(f"Skip {img_path}: {e}")

# Zielgröße
min_w = min(w for _, w, _ in sizes)
min_h = min(h for _, _, h in sizes)
print("Target size:", min_w, "x", min_h)

# Zuschneiden
for img_path, w, h in sizes:
    with Image.open(img_path) as img:
        left   = (w - min_w) // 2
        top    = (h - min_h) // 2
        right  = left + min_w
        bottom = top + min_h
        cropped = img.crop((left, top, right, bottom))
        cropped.save(out_dir / img_path.name)

print("Finished, all images cropped and saved in:", out_dir)
