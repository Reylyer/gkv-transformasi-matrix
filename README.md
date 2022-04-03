# gkv-transformasi-matrix
Tugas Kelas GKV Informatika Undip 2022

# Hasil
![muter](./media/videos/main/1080p60/TranslasiScene.gif)
![muter](./media/videos/main/1080p60/SkalaScene.gif)
![muter](./media/videos/main/1080p60/RotasiScene.gif)

# How to start?

## Setup Manim CE
Pertama kita perlu menginstall [manim community edition](https://github.com/ManimCommunity/manim). Ikuti cara install [disini](https://docs.manim.community/en/stable/installation.html).

Saya menyarankan menggunakan [chocolatey](https://chocolatey.org/install) karena mirip package manager di linux.
```
(admin) choco install manimce
(admin) choco install ffmpeg
```
## PASTI UDAH PUNYA PYTHON3 KAN!!!????
Kalo belom ni choco
```
(admin) choco install python --version=3.9.0
```
Ga disaranain ya ges make yang 3.10, banyak yang break mungkin manim juga termasuk

## Ngerender
Untuk merender pastiin udah di root repo ini lalu run
```
manim main.py
```
itu untuk full resolution ya ges bakal luama batttt. Ni yang kenceng no root low res
```
manim -ql main.py
```
```
-q, --quality [l|m|h|p|k]       Render quality at the follow resolution
                                  framerates, respectively: 854x480 15FPS,
                                  1280x720 30FPS, 1920x1080 60FPS, 2560x1440
                                  60FPS, 3840x2160 60FPS
```
