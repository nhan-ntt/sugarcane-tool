from sqlalchemy import Column, Integer, String, Text, ForeignKey, Index
from sqlalchemy.orm import relationship
from database import Base


class Genome(Base):
    __tablename__ = "genomes"

    # ID ngắn gọn: "R570", "AP85_441", "SSpon"
    id = Column(String, primary_key=True, index=True)
    name = Column(String)  # VD: Saccharum hybrid R570
    fasta_path = Column(String)  # Đường dẫn file: data/R570.fasta
    cds_path = Column(String)  # File CDS (.cds.fna) - MỚI
    protein_path = Column(String)  # File Protein (.faa) - MỚI

    # Quan hệ 1-N với Gene
    genes = relationship("Gene", back_populates="genome_ref")


class Gene(Base):
    __tablename__ = "genes"

    id = Column(Integer, primary_key=True, index=True)
    gene_id = Column(String, index=True)  # Không unique global nữa, chỉ unique trong 1 genome

    # --- THÊM CỘT NÀY ---
    genome_id = Column(String, ForeignKey("genomes.id"), index=True)
    genome_ref = relationship("Genome", back_populates="genes")
    # --------------------

    chromosome = Column(String, index=True)
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String)
    description = Column(Text)

    # Index tìm kiếm nhanh theo bộ gen và vị trí
    __table_args__ = (
        Index('idx_genome_loc', 'genome_id', 'chromosome', 'start', 'end'),
    )