from sqlalchemy import Column, Integer, String, Text, Index
from database import Base

class Gene(Base):
    __tablename__ = "genes"
    id = Column(Integer, primary_key=True, index=True)
    gene_id = Column(String, unique=True, index=True)
    chromosome = Column(String, index=True)
    start = Column(Integer, index=True)
    end = Column(Integer, index=True)
    strand = Column(String)
    description = Column(Text)
    __table_args__ = (Index('idx_region', 'chromosome', 'start', 'end'),)