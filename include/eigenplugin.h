template<typename OtherDerived>
inline Scalar frobInner(const MaxtrixBase<OtherDerived& other) const {
  return (derived().cwiseProduct(other.derived())).sum();
}
