from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
from typing import List

class AcademicFigureFormatter:
    """Padroniza figuras e cria painéis compostos para publicação."""

    def __init__(self, dpi=300, font_family="Arial"):
        self.dpi = dpi
        self.font_family = font_family

    def converter_para_alta_res(self, input_path: Path, output_path: Path):
        """Converte imagem para alta resolução (TIFF)."""
        img = Image.open(input_path)
        # Se for PNG com transparência, converte para RGB para TIFF
        if img.mode in ("RGBA", "P"):
            img = img.convert("RGB")
        
        output_path = output_path.with_suffix('.tiff')
        img.save(output_path, format='TIFF', dpi=(self.dpi, self.dpi), compression='tiff_lzw')
        return output_path

    def criar_painel(self, image_paths: List[Path], output_path: Path, labels: List[str], layout: str):
        """Cria um painel composto (ex: 1x2) com labels (A, B)."""
        rows, cols = map(int, layout.split('x'))
        
        images = [Image.open(p) for p in image_paths]
        
        # Assume que todas as imagens têm o mesmo tamanho para simplificar
        w, h = images[0].size
        
        panel_w = w * cols
        panel_h = h * rows
        
        panel = Image.new('RGB', (panel_w, panel_h), (255, 255, 255))
        draw = ImageDraw.Draw(panel)
        
        # Tentar carregar uma fonte, senão usa a padrão
        try:
            font = ImageFont.truetype(self.font_family, size=int(h * 0.05))
        except:
            font = ImageFont.load_default()

        for idx, img in enumerate(images):
            r = idx // cols
            c = idx % cols
            
            x_pos = c * w
            y_pos = r * h
            
            # Colar imagem
            if img.mode in ("RGBA", "P"):
                img = img.convert("RGB")
            panel.paste(img, (x_pos, y_pos))
            
            # Desenhar Label (A, B, C...)
            label = labels[idx]
            draw.text((x_pos + 10, y_pos + 10), label, fill=(0, 0, 0), font=font)

        output_path = output_path.with_suffix('.tiff')
        panel.save(output_path, format='TIFF', dpi=(self.dpi, self.dpi), compression='tiff_lzw')
        return output_path
