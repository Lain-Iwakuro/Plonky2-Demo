use anyhow::Result;
use plonky2::field::types::Field;
use plonky2::iop::target::Target;
use plonky2::iop::witness::{PartialWitness, WitnessWrite};
use plonky2::plonk::circuit_builder::CircuitBuilder;
use plonky2::plonk::circuit_data::CircuitConfig;
use plonky2::plonk::config::{GenericConfig, PoseidonGoldilocksConfig};
// use rand::random;
use rand_chacha::ChaChaRng;
use rand_chacha::rand_core::SeedableRng;
use rand::Rng;


/// An example of using Plonky2 to prove a statement of the form
/// "I know A * B".
fn main() -> Result<()> {
    const D: usize = 2;
    type C = PoseidonGoldilocksConfig;
    type F = <C as GenericConfig<D>>::F;

    let config = CircuitConfig::standard_recursion_config();
    let mut builder = CircuitBuilder::<F, D>::new(config);

    // The arithmetic circuit.
    let m = 64;
    let mut a = Vec::<Vec<Target>>::new();
    let mut b = Vec::<Vec<Target>>::new();
    let mut c = Vec::<Vec<Target>>::new();
    for i in 0..m {
        a.push(Vec::<Target>::new());
        b.push(Vec::<Target>::new());
        c.push(Vec::<Target>::new());
        (0..m).for_each(|_j: usize| {
            a[i].push(builder.add_virtual_target());
            b[i].push(builder.add_virtual_target());
        });
    }
    for i in 0..m {
        for j in 0..m {
            let initial = builder.constant(F::from_canonical_u32(0));
            let mut current_target = initial;
            for k in 0..m {
                let c_i_j_k = builder.mul(a[i][k], b[k][j]);
                current_target = builder.add(current_target, c_i_j_k);
            }
            c[i].push(current_target);
        }
    }

    // Public inputs are the initial value (provided below) and the result (which is generated).
    for i in 0..m {
        for j in 0..m {
            builder.register_public_input(a[i][j]);
            builder.register_public_input(b[i][j]);
            builder.register_public_input(c[i][j]);
        }
    }
    let mut pw = PartialWitness::new();
    let mut rand = ChaChaRng::from_entropy();
    for i in 0..m {
        for j in 0..m {
            pw.set_target(a[i][j], F::from_canonical_u32(rand.gen_range(u32::MIN..u32::MAX)));
            pw.set_target(b[i][j], F::from_canonical_u32(rand.gen_range(u32::MIN..u32::MAX)));
        }
    }
    let data = builder.build::<C>();
    let proof = data.prove(pw)?;

    println!("length of proof.public_inputs is {}", proof.public_inputs.len());
    /*
    let name = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    for k in 0..3 {
        println!("matrix {} = ", name[k]);
        for i in 0..m {
            for j in 0..m {
                print!("{:>3} ", proof.public_inputs[3 * (i * m + j) + k]);
            }
            println!("");
        }
    }
    println!("C should equals A * B");
    */
    data.verify(proof)
}
