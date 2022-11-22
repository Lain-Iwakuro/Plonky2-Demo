use alloc::vec;
use alloc::vec::Vec;
use core::marker::PhantomData;

use crate::field::extension::{Extendable, FieldExtension};
use crate::hash::hash_types::{HashOut, HashOutTarget, MerkleCapTarget, RichField};
use crate::hash::hashing::{HashConfig, PlonkyPermutation};
use crate::hash::merkle_tree::MerkleCap;
use crate::iop::ext_target::ExtensionTarget;
use crate::iop::target::Target;
use crate::plonk::circuit_builder::CircuitBuilder;
use crate::plonk::config::{AlgebraicHasher, GenericHashOut, Hasher};

/// Observes prover messages, and generates challenges by hashing the transcript, a la Fiat-Shamir.
#[derive(Clone)]
pub struct Challenger<F: RichField, HC: HashConfig, H: Hasher<F, HC>>
where
    [(); HC::WIDTH]:,
{
    pub(crate) sponge_state: [F; HC::WIDTH],
    pub(crate) input_buffer: Vec<F>,
    output_buffer: Vec<F>,
    _phantom: PhantomData<H>,
}

/// Observes prover messages, and generates verifier challenges based on the transcript.
///
/// The implementation is roughly based on a duplex sponge with a Rescue permutation. Note that in
/// each round, our sponge can absorb an arbitrary number of prover messages and generate an
/// arbitrary number of verifier challenges. This might appear to diverge from the duplex sponge
/// design, but it can be viewed as a duplex sponge whose inputs are sometimes zero (when we perform
/// multiple squeezes) and whose outputs are sometimes ignored (when we perform multiple
/// absorptions). Thus the security properties of a duplex sponge still apply to our design.
impl<F: RichField, HC: HashConfig, H: Hasher<F, HC>> Challenger<F, HC, H>
where
    [(); HC::WIDTH]:,
{
    pub fn new() -> Challenger<F, HC, H>
    where
        [(); HC::WIDTH]:,
    {
        Challenger {
            sponge_state: [F::ZERO; HC::WIDTH],
            input_buffer: Vec::with_capacity(HC::RATE),
            output_buffer: Vec::with_capacity(HC::RATE),
            _phantom: Default::default(),
        }
    }

    pub fn observe_element(&mut self, element: F)
    where
        [(); HC::WIDTH]:,
    {
        // Any buffered outputs are now invalid, since they wouldn't reflect this input.
        self.output_buffer.clear();

        self.input_buffer.push(element);

        if self.input_buffer.len() == HC::RATE {
            self.duplexing();
        }
    }

    pub fn observe_extension_element<const D: usize>(&mut self, element: &F::Extension)
    where
        F: RichField + Extendable<D>,
        [(); HC::WIDTH]:,
    {
        self.observe_elements(&element.to_basefield_array());
    }

    pub fn observe_elements(&mut self, elements: &[F])
    where
        [(); HC::WIDTH]:,
    {
        for &element in elements {
            self.observe_element(element);
        }
    }

    pub fn observe_extension_elements<const D: usize>(&mut self, elements: &[F::Extension])
    where
        F: RichField + Extendable<D>,
        [(); HC::WIDTH]:,
    {
        for element in elements {
            self.observe_extension_element(element);
        }
    }

    pub fn observe_hash<OHC: HashConfig, OH: Hasher<F, OHC>>(&mut self, hash: OH::Hash)
    where
        [(); OHC::WIDTH]:,
    {
        self.observe_elements(&hash.to_vec())
    }

    pub fn observe_cap<OHC: HashConfig, OH: Hasher<F, OHC>>(&mut self, cap: &MerkleCap<F, OHC, OH>)
    where
        [(); OHC::WIDTH]:,
    {
        for &hash in &cap.0 {
            self.observe_hash::<OHC, OH>(hash);
        }
    }

    pub fn get_challenge(&mut self) -> F
    where
        [(); HC::WIDTH]:,
    {
        // If we have buffered inputs, we must perform a duplexing so that the challenge will
        // reflect them. Or if we've run out of outputs, we must perform a duplexing to get more.
        if !self.input_buffer.is_empty() || self.output_buffer.is_empty() {
            self.duplexing();
        }

        self.output_buffer
            .pop()
            .expect("Output buffer should be non-empty")
    }

    pub fn get_n_challenges(&mut self, n: usize) -> Vec<F>
    where
        [(); HC::WIDTH]:,
    {
        (0..n).map(|_| self.get_challenge()).collect()
    }

    pub fn get_hash(&mut self) -> HashOut<F>
    where
        [(); HC::WIDTH]:,
    {
        HashOut {
            elements: [
                self.get_challenge(),
                self.get_challenge(),
                self.get_challenge(),
                self.get_challenge(),
            ],
        }
    }

    pub fn get_extension_challenge<const D: usize>(&mut self) -> F::Extension
    where
        F: RichField + Extendable<D>,
        [(); HC::WIDTH]:,
    {
        let mut arr = [F::ZERO; D];
        arr.copy_from_slice(&self.get_n_challenges(D));
        F::Extension::from_basefield_array(arr)
    }

    pub fn get_n_extension_challenges<const D: usize>(&mut self, n: usize) -> Vec<F::Extension>
    where
        F: RichField + Extendable<D>,
        [(); HC::WIDTH]:,
    {
        (0..n)
            .map(|_| self.get_extension_challenge::<D>())
            .collect()
    }

    /// Absorb any buffered inputs. After calling this, the input buffer will be empty, and the
    /// output buffer will be full.
    fn duplexing(&mut self)
    where
        [(); HC::WIDTH]:,
    {
        assert!(self.input_buffer.len() <= HC::RATE);

        // Overwrite the first r elements with the inputs. This differs from a standard sponge,
        // where we would xor or add in the inputs. This is a well-known variant, though,
        // sometimes called "overwrite mode".
        for (i, input) in self.input_buffer.drain(..).enumerate() {
            self.sponge_state[i] = input;
        }

        // Apply the permutation.
        self.sponge_state = H::Permutation::permute(self.sponge_state);

        self.output_buffer.clear();
        self.output_buffer
            .extend_from_slice(&self.sponge_state[0..HC::RATE]);
    }

    pub fn compact(&mut self) -> [F; HC::WIDTH] {
        if !self.input_buffer.is_empty() {
            self.duplexing();
        }
        self.output_buffer.clear();
        self.sponge_state
    }
}

impl<F: RichField, HC: HashConfig, H: AlgebraicHasher<F, HC>> Default for Challenger<F, HC, H>
where
    [(); HC::WIDTH]:,
{
    fn default() -> Self
    where
        [(); HC::WIDTH]:,
    {
        Self::new()
    }
}

/// A recursive version of `Challenger`. The main difference is that `RecursiveChallenger`'s input
/// buffer can grow beyond `HC::RATE`. This is so that `observe_element` etc do not need access
/// to the `CircuitBuilder`.
pub struct RecursiveChallenger<
    F: RichField + Extendable<D>,
    HC: HashConfig,
    H: AlgebraicHasher<F, HC>,
    const D: usize,
> where
    [(); HC::WIDTH]:,
{
    sponge_state: [Target; HC::WIDTH],
    input_buffer: Vec<Target>,
    output_buffer: Vec<Target>,
    __: PhantomData<(F, H)>,
}

impl<F: RichField + Extendable<D>, HC: HashConfig, H: AlgebraicHasher<F, HC>, const D: usize>
    RecursiveChallenger<F, HC, H, D>
where
    [(); HC::WIDTH]:,
{
    pub fn new(builder: &mut CircuitBuilder<F, D>) -> Self
    where
        [(); HC::WIDTH]:,
    {
        let zero = builder.zero();
        Self {
            sponge_state: [zero; HC::WIDTH],
            input_buffer: Vec::new(),
            output_buffer: Vec::new(),
            __: PhantomData,
        }
    }

    pub fn from_state(sponge_state: [Target; HC::WIDTH]) -> Self {
        Self {
            sponge_state,
            input_buffer: vec![],
            output_buffer: vec![],
            __: PhantomData,
        }
    }

    pub(crate) fn observe_element(&mut self, target: Target)
    where
        [(); HC::WIDTH]:,
    {
        // Any buffered outputs are now invalid, since they wouldn't reflect this input.
        self.output_buffer.clear();

        self.input_buffer.push(target);
    }

    pub fn observe_elements(&mut self, targets: &[Target])
    where
        [(); HC::WIDTH]:,
    {
        for &target in targets {
            self.observe_element(target);
        }
    }

    pub fn observe_hash(&mut self, hash: &HashOutTarget)
    where
        [(); HC::WIDTH]:,
    {
        self.observe_elements(&hash.elements)
    }

    pub fn observe_cap(&mut self, cap: &MerkleCapTarget)
    where
        [(); HC::WIDTH]:,
    {
        for hash in &cap.0 {
            self.observe_hash(hash)
        }
    }

    pub fn observe_extension_element(&mut self, element: ExtensionTarget<D>)
    where
        [(); HC::WIDTH]:,
    {
        self.observe_elements(&element.0);
    }

    pub fn observe_extension_elements(&mut self, elements: &[ExtensionTarget<D>])
    where
        [(); HC::WIDTH]:,
    {
        for &element in elements {
            self.observe_extension_element(element);
        }
    }

    pub fn get_challenge(&mut self, builder: &mut CircuitBuilder<F, D>) -> Target
    where
        [(); HC::WIDTH]:,
    {
        self.absorb_buffered_inputs(builder);

        if self.output_buffer.is_empty() {
            // Evaluate the permutation to produce `r` new outputs.
            self.sponge_state = builder.permute::<HC, H>(self.sponge_state);
            self.output_buffer = self.sponge_state[0..HC::RATE].to_vec();
        }

        self.output_buffer
            .pop()
            .expect("Output buffer should be non-empty")
    }

    pub fn get_n_challenges(&mut self, builder: &mut CircuitBuilder<F, D>, n: usize) -> Vec<Target>
    where
        [(); HC::WIDTH]:,
    {
        (0..n).map(|_| self.get_challenge(builder)).collect()
    }

    pub fn get_hash(&mut self, builder: &mut CircuitBuilder<F, D>) -> HashOutTarget
    where
        [(); HC::WIDTH]:,
    {
        HashOutTarget {
            elements: [
                self.get_challenge(builder),
                self.get_challenge(builder),
                self.get_challenge(builder),
                self.get_challenge(builder),
            ],
        }
    }

    pub fn get_extension_challenge(
        &mut self,
        builder: &mut CircuitBuilder<F, D>,
    ) -> ExtensionTarget<D>
    where
        [(); HC::WIDTH]:,
    {
        self.get_n_challenges(builder, D).try_into().unwrap()
    }

    /// Absorb any buffered inputs. After calling this, the input buffer will be empty, and the
    /// output buffer will be full.
    fn absorb_buffered_inputs(&mut self, builder: &mut CircuitBuilder<F, D>)
    where
        [(); HC::WIDTH]:,
    {
        if self.input_buffer.is_empty() {
            return;
        }

        for input_chunk in self.input_buffer.chunks(HC::RATE) {
            // Overwrite the first r elements with the inputs. This differs from a standard sponge,
            // where we would xor or add in the inputs. This is a well-known variant, though,
            // sometimes called "overwrite mode".
            for (i, &input) in input_chunk.iter().enumerate() {
                self.sponge_state[i] = input;
            }

            // Apply the permutation.
            self.sponge_state = builder.permute::<HC, H>(self.sponge_state);
        }

        self.output_buffer = self.sponge_state[0..HC::RATE].to_vec();

        self.input_buffer.clear();
    }

    pub fn compact(&mut self, builder: &mut CircuitBuilder<F, D>) -> [Target; HC::WIDTH] {
        self.absorb_buffered_inputs(builder);
        self.output_buffer.clear();
        self.sponge_state
    }
}

#[cfg(test)]
mod tests {
    use crate::field::types::Sample;
    use crate::iop::challenger::{Challenger, RecursiveChallenger};
    use crate::iop::generator::generate_partial_witness;
    use crate::iop::target::Target;
    use crate::iop::witness::{PartialWitness, Witness};
    use crate::plonk::circuit_builder::CircuitBuilder;
    use crate::plonk::circuit_data::CircuitConfig;
    use crate::plonk::config::{GenericConfig, PoseidonGoldilocksConfig, PoseidonHashConfig};

    #[test]
    fn no_duplicate_challenges() {
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type HCO = PoseidonHashConfig;
        type HCI = HCO;
        type F = <C as GenericConfig<HCO, HCI, D>>::F;
        let mut challenger =
            Challenger::<F, HCI, <C as GenericConfig<HCO, HCI, D>>::InnerHasher>::new();
        let mut challenges = Vec::new();

        for i in 1..10 {
            challenges.extend(challenger.get_n_challenges(i));
            challenger.observe_element(F::rand());
        }

        let dedup_challenges = {
            let mut dedup = challenges.clone();
            dedup.dedup();
            dedup
        };
        assert_eq!(dedup_challenges, challenges);
    }

    /// Tests for consistency between `Challenger` and `RecursiveChallenger`.
    #[test]
    fn test_consistency() {
        const D: usize = 2;
        type C = PoseidonGoldilocksConfig;
        type HCO = PoseidonHashConfig;
        type HCI = HCO;
        type F = <C as GenericConfig<HCO, HCI, D>>::F;

        // These are mostly arbitrary, but we want to test some rounds with enough inputs/outputs to
        // trigger multiple absorptions/squeezes.
        let num_inputs_per_round = vec![2, 5, 3];
        let num_outputs_per_round = vec![1, 2, 4];

        // Generate random input messages.
        let inputs_per_round: Vec<Vec<F>> = num_inputs_per_round
            .iter()
            .map(|&n| F::rand_vec(n))
            .collect();

        let mut challenger =
            Challenger::<F, HCI, <C as GenericConfig<HCO, HCI, D>>::InnerHasher>::new();
        let mut outputs_per_round: Vec<Vec<F>> = Vec::new();
        for (r, inputs) in inputs_per_round.iter().enumerate() {
            challenger.observe_elements(inputs);
            outputs_per_round.push(challenger.get_n_challenges(num_outputs_per_round[r]));
        }

        let config = CircuitConfig::standard_recursion_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let mut recursive_challenger =
            RecursiveChallenger::<F, HCI, <C as GenericConfig<HCO, HCI, D>>::InnerHasher, D>::new(
                &mut builder,
            );
        let mut recursive_outputs_per_round: Vec<Vec<Target>> = Vec::new();
        for (r, inputs) in inputs_per_round.iter().enumerate() {
            recursive_challenger.observe_elements(&builder.constants(inputs));
            recursive_outputs_per_round.push(
                recursive_challenger.get_n_challenges(&mut builder, num_outputs_per_round[r]),
            );
        }
        let circuit = builder.build::<HCO, HCI, C>();
        let inputs = PartialWitness::new();
        let witness = generate_partial_witness(inputs, &circuit.prover_only, &circuit.common);
        let recursive_output_values_per_round: Vec<Vec<F>> = recursive_outputs_per_round
            .iter()
            .map(|outputs| witness.get_targets(outputs))
            .collect();

        assert_eq!(outputs_per_round, recursive_output_values_per_round);
    }
}
