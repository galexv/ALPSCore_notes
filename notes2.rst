Notes on ALPSCore development
#############################

Custom types in Accumulators
============================

The ``mean_type<A>`` "meta-returns" ``double`` if ``value_type<A>``
returns an integral type (that is,
``boost::is_integral<value_type<A>::type`` is a true type); otherwise,
it is defined to return ``value_type<A>``. So, for an accumulator of
``custom_type`` the ``mean_type`` returns the ``custom_type`. Ok.

Next: tries to convert ``B::count()`` result to
``alps::hdf5::scalar_type<mean_type>>::type``. So, we need a
specialization of the trait ``scalar_type<T>``, which is in
``alps/hdf5`` include directory --- used for serialization. The
implicit assumption here is that ``scalar_type<T>`` (let's call it
``ST``) is something that can hold the count, and that can be a
divisor for type T to obtain the mean, which is of type T itself. That
is, ``ST x=int()`` is valid, and ``T==T/ST``. 

Next: ``operator<<`` must be defined, to output the value of
``custom_type``.

Next: ``custom_type operator*(const custom_type&, const ST&)`` should
be defined, with the semantics of ``custom_type sum =
custom_type(mean) * ST`` (that is, multiplication should be inverse of
division).

Suggestion: instead of ``hdf5::scalar_type`` (which has a meaning of
"underlying data storage type") there should be
``accumulators::scalar_type`` or ``numeric::scalar_type`` with the
meaning of "type that can be used for invertible scaling". The
corresponding test would be the invertibility of ``*`` and ``/``.

Next: serialization (for MPI). At minimum, we need to iterate over all
elements of the ``custom_type``.  (Otherwise, we should specialize the
trait ``hdf5::is_continuous<???>`` --- and other traits).

How MPI (in this case, ``collective_merge()``) works? Suppose, the
data type is ``T`` (``custom_type`` in my case). Suppose, ``typedef
hdf5::scalar_type<T>::type ST``.

    # ``reduce_if(comm, arg=T(value), res=value, op=std::plus<ST>(), root)`` is
    called.
    # The ``ST`` type should be C++ scalar (``boost::is_scalar<ST>``
    is true_type), otherwise exception is thrown. 
    # ``alps::mpi::reduce(comm, arg, res, op, root`` is called.
    # ``reduce_impl(comm,in_val=arg,out_val=res,op,root, 
    is_scalar=boost::is_scalar<T>::type(),
    is_cont=hdf5::is_content_continuous<T>::type() )`` is called.
    # If ``T`` is scalar, then ``is_cont`` is ignored and
    ``boost::mpi::reduce()`` is called.
    # If ``is_scalar==false_type``, but ``is_cont==true_type`` (continuous
    type), then:
    
        # ``hdf5::get_extent(in_val)`` is used to inquire the extents;
        # ``hdf5::set_extent(out_val, ext)`` is used to set the target
        extents;
        # ``boost::mpi::get_mpi_datatype(ST())`` is used to inquire
        the MPI equivalent of the underlying scalar datatype;
        # ``hdf5::get_pointer(in_val)``,
        ``hdf5::get_pointer(out_val)`` are used to get the pointers to
        the continuous data areas;
        # ``MPI_Reduce()`` is called.

    # If ``is_scalar==false_type`` and ``is_cont==false_type``
    (discontinuous type), then:
    
        # ``hdf5::is_vectorizable(in_val)`` is inquired. If false, an
        exception is thrown;
        # Like above, ``get_extent(in_val)`` etc.
        # Allocate ``std::vector<ST>`` of the proper length.
        # call ``copy_to_buffer(....)``
        # Like above, call ``set_extent(out_val,...)``.
        # call ``copy_from_buffer(...)``.

Those ``copy_to_buffer()``, ``copy_from_buffer()`` work recursively on
vectorizable objects down to vectors of scalars.

In our particular case, Eigen's matrices has continuous storage --- so,
we should just declare ``custom_type`` to be continuous (otherwise, it
is expected to be vectorizable).

Next: implement ``save()`` and ``load()``.

Next: Arithmetic operators. It is possible to declare them as friends
to a base class, define them there (to throw an exception), then
inherit from this class. It works, but then it would probably be not
possible to redefine the operators to do anything sensible.

``mean.hpp`` defines arithmetic operators with scalar RHS by casting
the scalar to ``alps::element_type<mean_type>::type``. It is
confusing: it assumes that the data type in question consists of
elements of the math-scalar type --- which is not necessarily the
case. It should be, again, ``numeric::scalar_type``. Generally, all
uses of ``element_type`` are suspects for replacing to
``numeric::scalar_type``!

The ``element_type<T>`` metafunction returns ``T::value_type`` if
defined, otherwise ``T``. For now, let's just introduce
``value_type`` in the ``custom_type``.

Next: we have a call like ``base_wrapper<custom_type>(lhs) AUGOP
base_wrapper<custom_type::value_type>(rhs)``, which is not defined ---
because the metafunction ``wrap_value_type<X>`` does not return the
correct result for our ``custom_type``. We have to redefine the
metafunction to return ``base_wrapper<numeric::scalar_type<X>::type>``
(for now, let it use ``element_type`` instead). This can be achieved
by checking for non-scalar type (that is, the type ``T`` whose
``scalar_type<T>`` is not the same as ``T``).

Problem: ``alps::has_value_type<custom_type>`` returns false_type,
because, apparently, when it is first instantiated, ``custom_type`` is
declared but not yet defined, and therefore does NOT have
``value_type`` member. Workaround: the custom type should be fully
defined before including any ALPSCore files. Inconvenient. :(

We might be able to circumvent this by postponing the corresponding
instantiations.  

Next: Problems with ``Result<...>::scalar_result_type``. We should do
the same as we did with ``wrap_value_type``: let it return the type
based on ``numeric::scalar_type`` (which is ``element_type`` for
now). But i will run into problems with ``vector_result_type``: who
said that "non-scalar" is a vector?? (This is fixed now).

Sometimes we need operators (for ``custom_type x; scalar_type s;``)
like ``s+x`` and ``s-x``. Let us assume and require that the custom
type is mathematically sane, that is:

   s + x =  x + s
   s - x = -x + s
   s * x = x * s
   s * x * y = x * y * s

Question: how at all error bar behaves when multiplying
non-commutative values?

Next: ``alps::numeric::inf<custom_type>``: the notion of the "infinite
value". This template must be specialized in the ``alps::numeric``
namespace. It is returned when the error bar cannot be determined
because there are not enough bins. 

It is also returned when there is not enough bins to establish
the autocorrelation value --- which brings up a question: is
autocorrelation value for a matrix-valued function also a matrix?

The ``alps::numeric::inf<custom_type>`` is expected to be
default-constructible, and implicitly-convertible to ``custom_type``.

Next: ``alps::numeric::set_negative_0(T&)``. ALPSCore mistakes
``custom_type`` for a sequence (``is_sequence<custom_type>`` is
true_type), because of ``custom_type::value_type``. It has to be
explicitly declared as not a sequence.

This ``set_negative_0()`` sets its argument to 0 if it is
negative. What does it mean for a matrix to be negative, and why is it
needed? Answer: for matrices, we are doing element-by-element
operations.

Problem: pointers to the math functions, like ``sin(custom_type)``, if
declared as friend-functions to a parent class, cannot be looked up
properly. Therefore, we have just to move them in the namespace
proper. Moreover, those functions must expect ``custom_type`` rather
than ``const custom_type&`` --- which is OK, since most probably a
copy will be needed... (Needs checking: is it optimized properly?)

Next: bumped into the problem that the jacknife machinery is private. 
Inside the call of
``Result<custom_type,...>::transform(...Result<double,...> arg...)``
it tries to call ``arg.generate_jacknife()``, which is private. Let's
just make it public for now, and hide member variables behind
accessors. 

There are numerous places, related to MPI reduction, where
``element_type<custom_type>`` or ``hdf5::scalar_type<custom_type>`` is
used. It is a tricky question which needs further investigation: can
it always be replaced with ``numeric::scalar<custom_type>``? In
practical terms, for now these two types must (most probably) be the
same, because, the way MPI reduction is implemented, sum of two
objects is a sum of its elements. On the other hand,
``numeric::scalar<T>`` is used for scaling (multiplication &
division), while ``hdf5::scalar_type<T>`` is involved in MPI reduction
(summation) --- so, they can be distinct.

Semantics (all are in namespace ``alps::``):

    * ``hdf5::scalar_type<T>``: the "constituent" type of T; the type of the
    elements of the underlying data structure. Used for saving/loading
    and MPI reduction. The same as ``T`` by default.

    * ``is_sequence<T>``: true_type if ``T`` is a sequence. By
      default, true if ``T::value_type`` exists.

    * ``numeric::scalar<T>``: the "mathematical scalar" of T; the type
    of the thing that can be used to scale values of T. Used for
    arithmetic operations on T. I made it the same as
    ``element_type<T>`` by default.

    * ``average_type<T>``: the type that can hold average value of
      several ``T`` values.

    * ``element_type<T>``: defined as ``T::value_type`` if it exists;
    otherwise, it's ``T``. Used to define slices of sequences and the
    covariance type. Access to a slice of a sequence returns its
    ``element_type``.

    * ``covariance_type<T>``: matrix of average types of
    ``element<T>::type`` if ``T`` is a sequence, otherwise
    ``average_type<T>``. Outer product of two sequences returns
    ``covariance_type``.
      
This means that if ``custom_type`` is not a sequence, its
``element_type`` must be ``custom_type``!


Multiplication ``operator*(custom_type, custom_type)`` is used for
error calculation (in ``error`` feature and in ``binning_analysis``
feature). How should it behave for matrix types? Emanuel's answer:
"always element-wise". In this case, should all supported type be
"sequences"?

If ``custom_type`` is a sequence, it is used to define
``set_negative_0()`` (as a loop over elements of the sequence), and
``checked_divide()`` (as a loop). Also, the following traits are to be
defined for the sequence type ``T``:

    * Function ``size(const T& seq)`` (default: ``seq.size()``);

    * Metafunction ``covariance_type<T>`` (default:
      ``boost::numeric::ublas::matrix`` of average-type of ``element_type<T>::type``);

    * Function ``covariance_type<T>::type outer_product(T a, T b)``
      (default: UBlas matrix, outer product of ``a`` and ``b``
      regarded as vectors); used in ``max_num_binning.hpp``, in
      ``Result<...>::covariance()`` and in
      ``Result<...>::accurate_covariance()``.

    * Metafunction ``slice_index<T>`` (default: ``std::size_t``);

    * Function ``std::pair<slice_index<T>::type,
      std::pair<slice_index<T>::type> slices(const T& seq)`` (default:
      returns a pair of 0, ``seq.size()``); used in
      ``binning_analysis.hpp`` in
      ``Accumulator<...>::converged_errors()`` method.

    * Function ``element_type<T>::type slice_value(const T& seq, unsigned i)``
      (default: checks the size, returns `seq[i]`` or default element
      value); used in ``binning_analysis.hpp`` in
      ``Accumulator<...>::converged_errors()`` method.

    * Function ``element_type<T>::type& slice_value(T& seq, unsigned i)``
      (default: returns ``seq[i]``); used in ``binning_analysis.hpp``
      in ``Accumulator<...>::converged_errors()`` method.

    * Functor class ``slice_it<T>`` with the method 
      ``element_type<T>::type operator(const T& seq, slice_index<T>::type i)``
      (default: returns ``slice_value(seq,i)``). Is not used anywhere.

    * Function ``T checked_divide(T a, const T& b)`` (default: element-wise
      ``checked_divide()``); not used anywhere.
    
    * Function ``void check_size(T& a, const T& b)`` (default: do
      ``a.resize(b.size())``, specialized for vectors to resize a
      zero-sized ``a``, otherwise throw); used in several places,
      e.g. when a value is added to an accumulator.

    * Function ``void set_negative_0(T&)`` (default: set each negative
      element to 0); used in ``binning_analysis.hpp``
      in ``Accumulator<...>::autocorrelation()`` method.

The default implementations of these traits assume that ``seq.size()``
returns the size of the sequence, and that ``seq[i]`` returns the i-th
element of the sequence. 

It seems that an Eigen's matrix should be a sequence, then: especially
because of a covariance type.

Interestingly, ``operator+=(custom_type,scalar)`` is not needed for
the program to compile. Apparently, the corresponding operation on
results is defined via ``operator+``. 


Relationships between various types
===================================

Named accumulators: ``FullBinningAccumulator<double>``.
    * Can be added to ``accumulator_set`` (aka ``impl::wrapper_set<accumulator_wrapper>``).
    * Contains a name
    * Contains ``wrapper`` holding ``shared_ptr<accumulator_wrapper>``.
    * Contains a type ``accumulator_type`` which is a type of the
      corresponding "raw accumulator".
    * Contains a type ``result_type`` which is a type of the
      corresponding "raw result".

The mapping between the "named" accumulators and "raw" accumulators'
feature tags::
  MeanAccumulator : mean_tag, count_tag
  NoBinningAccumulator : error_tag, mean_tag, count_tag
  LogBinningAccumulator : binning_analysis_tag, error_tag, mean_tag, count_tag
  FullBinningAccumulator : max_num_binning_tag, binning_analysis_tag, error_tag, mean_tag, count_tag
      
Accumulator|Result sets:
    * Basically, a map containing shared pointers to ``accumulator_wrapper`` | ``result_wrapper``.
    * Access operator ``operator[](name)`` returns a reference to ``accumulator_wrapper`` | ``result_wrapper``.

Accumulator|Result wrappers (``accumulator_wrapper`` | ``result_wrapper``):
    * Contains a variant of ``shared_ptr< base_wrapper<T> >``, where ``T`` runs over all supported data types,
    * The pointer actually points to ``derived_accumulator_wrapper<A>`` | ``derived_result_wrapper<A>``.
      `QUESTION:` where does it take its value from??
    * Supports ``mean()``, ``error()`` methods, as well as arithmetic methods.
    * The method calls are forwarded via virtual methods of ``base_wrapper<T>``
      to the object actually held in the variant.
    * Supports method ``A& extract<A>()`` returning underlying raw accumulator|result ``A``.
    * There are free functions ``extract<A>()``.

The ``base_wrapper<T>`` class (where ``T`` is one of the supported
data types) inherits from the chain of ``impl::BaseWrapper<T,F,B>``
where ``F`` is a "feature tag" and ``B`` is a base class (next in
chain). At the end of the chain is ``detail::value_wrapper<T>``, which
contains only definition of ``value_type`` as ``T``.

Both ``derived_accumulator_wrapper<A>`` and
``derived_result_wrapper<A>`` inherit from ``derived_wrapper<A>``.

The ``derived_wrapper<A>`` class inherits from the chain of
``impl::DerivedWrapper<T,F,B>`` where ``F`` is a "feature tag" and
``B`` is a base class (next in chain). At the end of the chain there
is ``detail::foundation_wrapper<A>``, which inherits from
``base_wrapper<value_type<A>::type>``.

To check whether a given value of ``accumulator_wrapper`` related to a
particular named accumulator::
    using namespace alps::accumulators;
    #define SomeNamedAccumulator FullBinningAccumulator<double> /* or some other named accumulator */
    // Suppose we have a set of measurements mset:
    accumulator_set mset;
    mset << SomeNamedAccumulator("name");
    // Obtain the wrapper object by name:
    const accumulator_wrapper& acc=mset["name"]; 
    // Try to extract the raw accumulator:
    try {
        acc.extract<SomeNamedAccumulator::accumulator_type>();
    } catch (std::bad_cast&) {
        std::cerr << "acc does not wrap SomeNamedAccumulator or its subclass";
    }

   
    

Adding a new method to ``result_wrapper``
=========================================

How does ``result_wrapper::error<U>()`` work?
---------------------------------------------

Methods like ``result_wrapper::error<U>()`` use visitor to operate on
value of type ``shared_ptr<X>`` determined by meta-predicate
``detail::is_valid_argument< error_type<X>::type, U>``, where ``X`` is
``base_wrapper<T>``.  The ``error_type<X>`` is defined (as the
"identical" subclass of ``mean_type<X>``) in ``error.hpp``. The
``mean_type<X>`` is ``value_type<X>::type`` (if the latter type is not
integral, otherwise ``double``).

The metafunction ``value_type<X>`` returns ``X::value_type``. 

Therefore, the metafunction ``error_type<X>`` returns
``X::value_type``, which is inherited through the chain of
``impl::BaseWrapper<T,F,B>`` from ``detail::value_type<T>`` to be
``T``.

The meta-predicate ``detail::is_valid_argument<error_type<X>::type, U>``
therefore compares between two datatypes: ``detail::is_valid_argument<T,U>``,
where ``T`` is the underlying datatype of the ``result_wrapper<T>``,
and ``U`` is the template argument in the call ``error<U>()``. 

How can we add ``result_wrapper::autocorrelation<U>()``?
--------------------------------------------------------

The visitor classes are generated using macro
``ALPS_ACCUMULATOR_PROPERTY_PROXY(method,return_type)``. In our case,
the ``return_type`` is ``autocorrelation_type``; this macro then
generates call to metafunction ``autocorrelation_type``
(``autocorrelation_type<X>::type``). Therefore, we need a metafunction
``autocorrelation_type<X>`` that returns ``X::autocorrelation_type``
-- and it is there, in ``binning_analysis.hpp``, and is visible. This part works.

Then, ``X::autocorrelation()`` method is called (that is,
``base_wrapper<T>::autocorrelation()``), which compiles, finding the
(pure virtual!)  method in one of its parent classes,
``impl::BaseWrapper<T,binning_analysis_tag,B>``. It must get
dispatched to ``derived_result_wrapper<A>::autocorrelation()``, which
finds the method in one of its parent classes,
``impl::DerivedWrapper<T,binning_analysis_tag,B>``.


The method
``impl::DerivedWrapper<T,binning_analysis_tag,B>::autocorrelation()``
calls ``detail::autocorrelation_impl(arg)``, where ``arg`` is taken
from its ``m_data`` field; the latter is a protected member of type
``A`` inherited from ``detail::foundation_wrapper<A>``. The free
templated function ``detail::autocorrelation_impl<A>(const A& acc)``
calls ``autocorrelation(acc)`` (or throws if meta-predicate
``has_feature<A, binning_analysis_tag>`` is false). The
``autocorrelation(acc)`` simply calls ``acc.autocorrelation()``.


Adding serialization (save/load) capabilities.
==============================================

The class has to implement several concepts (the info is taken from
[https://alps.comp-phys.org/trac/wiki/NGSHDF5#Non-IntrusiveSerialization]
and from looking at ALPSCore code). The following template specializations must be defined for the type T (in ``alps::hdf5`` namespace):

    * Metafunction ``scalar_type<T>`` 
    * Meta-predicate ``has_complex_elements<T>`` 
    * Meta-predicate ``is_continuous<T>`` 
    * Function ``void save<T>(
      alps::hdf5::archive& ar, std::string const& path, const T& value,
      std::vector<std::size_t> size = std::vector<size_t>(),
      std::vector<std::size_t> chunk = std::vector<size_t>(),
      std::vector<std::size_t> offset = std::vector<size_t>())``
    * Function ``void load<T>(
      alps::hdf5::archive& ar, std::string const& path, T & value,
      std::vector<std::size_t> chunk = std::vector<size_t>(),
      std::vector<std::size_t> offset = std::vector<size_t>())``

in ``alps::hdf5::detail`` namespace, template specializations with static ``apply()`` method:

    * ``std::vector<std::size_t> get_extent<T>::apply(const T& value)``
    * ``void set_extent<T>::apply(T& value, const std::vector<std::size_t>& extent)``
    * ``bool is_vectorizable<T>::apply(const T& value)``
    * ``some_pointer_type get_pointer<T>::apply(T& value)``
    * ``some_const_pointer_type get_pointer<const T>::apply(const T& value)``

This is a mess --- this is documented poorly.

``save()`` ultimately calls ``alps::hdf5::archive::write<U>(path,
const U* val_addr, size_vec, chunk_vec, offset_vec)`` which is defined
only for native HDF5 types ``U``. Therefore,
``get_pointer<custom_type>()`` should return a native type pointer.

True-valued ``is_continuous<T>`` apparently means that if we have this::
    T* array_of_values=new T[N];
    U* ptr=get_pointer<T>::apply(array_of_values);
then ``ptr`` points to a *continuous* area of data associated with
``array_of_values``.


What to do with it?
-------------------

We can introduce "intrusive" ``save()`` and ``load()``. It works, as
long as no name is assigned when saving (as in
``ar[""]<<value``). This is because, e.g., feature= ``mean_tag``
expects that the mean value is "data" (not "datagroup"!) and is of the
type ``hdf5::scalar_type<mean_type>::type``, and also either scalar
(according to ``boost::is_scalar``) or ``get_extent(T()).size()``
returns the number of dimensions associated with the HDF5 data (``T``
is the datatype held in the accumulator).

Another requirement is that
``accumulators::accumulator_set::register_serializable_type<raw_acc_type>()``
must be called (once! before loading!), where ``raw_acc_type`` is the
"raw accumulator type" holding ``T``.

It is still a mess: when loading, each of the registered accumulator
types are tried in turn, until the first one capable of reading this
data is found. That results in ``double`` accumulator thinking that it
can read my ``custom_type`` accumulator, the ``double`` -holding
accumulator being constructed, and a failure when the ``double``
accumulator is converted to ``custom_type``.

Proposed solution: each registered type should announce an attribute
value that will be saved with the data; the loader should match the
attribute with the corresponding accumulator reader. It looks like a
not so big extension of ``can_load()`` method.

So, how does it work? ``wrapper_set<W>``, where ``W`` is
``accumulator_wrapper`` or ``results_wrapper``, contains
``std::vector<boost::shared_ptr<detail::serializable_type<W> > >
m_types`` member: a vector of shared pointers to
``detail::serializable_type<W>``. The latter is an abstract class,
from which ``detail::serializable_type_impl<W,A>`` is derived, where
``A`` is a raw accumulator type (the template parameter of the
``accumulators::accumulator_set::register_serializable_type<A>()``
above). The abstract class defines the following interface:

    * ``std::size_t rank() const`` : returns the "rank" of the loader;
      higher ranks are tried first. The default implementation
      behavior is to call ``A::rank()``, which returns 1 for the
      count-accumulator, and ``B::rank()+1`` for all other raw
      accumulators.
    * ``bool can_load(hdf5::archive&) const`` : returns whether this
      type can be loaded. Default implementation: returns
      ``A::can_load()``.
    * ``W* create(hdf5::archive&) const`` : returns a pointer to a
      newly-created object on a heap. Default implementation: ``return
      new W(A())``; that is, creates a new wrapper initialized with
      default-constructed raw accumulator.

The ``register_serializable_type()`` function inserts the new
serializable type in the vector (``m_types`` above) according to the
rank, maintaining the vector sorted in descending order (higher ranks
first).

So, for a user-defined type ``T``, i can either expand the default
implementation, or substitute it with another, specialized for each
``A<T>``. How would the alternative implementation look?

To load an accumulator, i still need to have a simple data (vector or
scalar) written in HDF5. This implementation allows for customized
accumulator (type ``A``) reading, while i rather need a customize data
type (type ``T``) reading!

Let's introduce a ``bool can_load(hdf5::archive& ar, const
std::string& name, size_t dim)`` as a trait for
``T``. (Is there already something like that? Does not look like.)
This will not work, because ``tau/data`` is an array of dimension 1
more than the saved data: it cannot be (naturally) generalized to an
arbitrary data. OTOH, i can explicitly save/load a vector of
general values --- why not?

So, do i need to extend the archive concepts to introduce
``can_load()``? 

For now, let's just have ``can_load()`` that checks name, datatype,
dimensionality (0 means check whether the type ``T`` and the saved
type are both scalars or both not).

mean value: (the tested type T is actually mean_type).
  1. Is data.
  2. Contains the datatype of the hdf5-scalar type.
  3. If T is boost-scalar, then data is scalar
  4. If T is not boost-scalar, then data is not scalar and has
     dimensions of ``get_extent(T()).size()``.

So, the call should be:
``trait<T>::can_load(ar,name,boost::is_scalar<T>.value?0:get_extent(T()).size())``.

Now the trait is coded and used. It has to be specialized for
``custom_type``, and has to work for the array of ``custom_type`` too.

When i am trying to save a vector of ``custom_type``, the traits like
``get_extent<T>`` and ``get_pointer<T>`` are called, rather than
``save()``/``load()``. How to prevent this? How does ``std::vector``
saves its elements?

It is enough to declare the type as non-continuous and
non-vectorizable; then ``save()`` or ``load()`` will be called
element-wise. However, it would break MPI functionality!

A compromise variant: the type can be declared non-continuous and
non-vectorizable, but content continuous --- will not work: it breaks
when MPI is used for a vector of ``custom_type``.

I was wrong about semantics of ``is_continuous<T>``: it returns
``true`` if an *array* of T is continuous data. 

How does ``hdf5::archive::read(string path, T* val, vector<size_t> chunk,
vector<size_t> offset)`` work?

  1. ``vector<size_t> data_size=extent(path)``:  sets the
     ``data_size`` vector to the extent of the saved data, which is:
    * For non-existing data: ``{0}``
    * For scalars: ``{1}``
    * For N-dim non-scalars: ``{0, ..., 0}`` (of length N), which then
      apparently gets filled with the actual dimensions of the data or
      attribute ``path``. 
  2. Verify that the number of elements in ``data_size`` (that is, the
     dimensionality of the data) is the same as the number of elements
     in ``chunk`` and ``offset``.
  3. For each dimension ``i``, ``data_size[i]`` should be not less
     than ``chunk[i]+offset[i]``. 
  4. String type, then each native type, are tried in order:
     1. Length is a multiple of ``chunk`` array elements.
     2. Raw array of the corresponding length is allocated.
     3. If chunks (that is, ``chunk`` array elements) are all the same
        as data dimensions (that is, ``data_size`` array elements),
        the data is read into the raw array;
     4. Otherwise:
        1. A hyperslab is selected, with each dimension ``i``
           starting at ``offset[i]`` and having ``chunk[i]`` elements.
        2. A "file dataspace" is created based on the hyperslab.
        3. A "memory dataspace" is created of dimensions from ``chunk`` array.
        4. The file dataspace is read into raw array described by the
           memory dataspace.
        5. `In other words`, a slab (a multi-dimensional sub-array) is
           cut out from the file data; the dimensions (sizes) of the slab are
           given by ``chunk[]``, and the offsets (positions) in the
           data structure in the file are given by ``offset``.
     5. The raw array is cast to the return value.

By implementing ``get_extent`` and ``save()``/``load()`` properly we
made the saving/loading of the custom type and vectors of custom type
work. 

A saved value or vector of custom type looks exactly like any other
scalar or vector in HDF5; to be precise, ``my_custom_type<T>`` looks
like a one-dimensional one-element vector of ``T``. How to distinguish
the custom type values?  For example, setting an attribute
``c++_type`` with value ``my_custom_type``. A loader can verify the
attribute and the type of the data entity (using ``bool
hdf5::archive::is_datatype<T>(string path)``) and decide whether the
data is loadable.

Because it appears to be perfectly fine to read my custom type into a
vector of any type, load attempts for accumulators should begin with
the custom type and only then try regular types.






Tests for proper accumulator statistics.
========================================

The errorbars (and correlation lengths) returned by binning
accumulators depends on the data generators, and differ widely from
non-binned estimators in cases of non-uniformly-distributed data. How
to (easily) introduce accumulator-specific errorbar verifiers?

One option, little bit verbose, is to ask to supply a predicate
functor: it will be given the errorbar as argument, and should return
true if the errorbar is acceptable, false otherwise.

Another option is to subclass the test class, declare a separate test
set, and code an errorbar tester explicitly --- and this is what i did.

Accumulator printing.
=====================

I want something like this: ``cout << fullprint(acc)``. That is,
``fullprint(acc)`` should return something which is printable: e.g.,
``detail::accumulator_print_proxy``. Possible, but impractical for
now.


Merging accumulators.
=====================

To make it finally work, i need "universal summation": vectors,vectors
of vectors... etc. It is possible (i made a sample program). However,
we need to figure out what defined in ALPSCore in terms of what.

``boost_array_functions.hpp``: 
Defines ``numeric::operator+=`` and the like for ``boost::array``. Not
relevant.

``functional.hpp``: defines classes ``numeric::plus<T,U,R>`` and
``numeric::plus<T,T,T>`` and the like in terms of
``boost::numeric::operators::operator+`` (and the like), which
presumably do element-wise vector operations. Some of them are using
``alps::numeric::operator+`` instead! (commit 66ba01bd).

``vector_functions.hpp``: defines ``operator+=`` and the like,
originally in terms of ``std::plus`` and the like. We should rewrite
them in terms of ``numeric::plus``.

Problem: accumulators have different number of bins. We need a "merge"
operation on vectors, with resizing the smaller vector to the larger
one. How would we do it?

Basically, it is "OP-EQ" operation.
If the left vector is smaller:::
  lsz=left.size();
  rsz=right.size();
  if (lsz<rsz) left.resize(rsz);
  left += right;

However, if the right vector is smaller, the operation becomes more
involved. We need either to copy (and then discard) the right vector
before resizing, or just introduce a special "merge":::
  lsz=left.size();
  rsz=right.size();
  if (lsz<rsz) left.resize(rsz); // now left is at least as big as right
  std::transform(right.begin(), right.end(), 
                 left.begin(),
                 left.begin(),
                 alps::numeric::plus<T,T,T>());

Looks like working.

Next problem: Half-finished merge of ``max_num_binning`` accumulator,
which does nothing except calling the merge of its base class (which
is ``binning_analysis`` accumulator), passes the test. We need a better
merge test!


