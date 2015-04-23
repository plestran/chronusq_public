template<typename OtherDerived>
inline Scalar frobInner(const MatrixBase<OtherDerived>& other) const {
  return (derived().cwiseProduct(other.derived())).sum();
}
