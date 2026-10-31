# ==============================================================================
# Estágio 1: Compilação Nativa do Motor Astronômico em Rust (`aa56-cli`)
# ==============================================================================
FROM rust:1.83-slim AS rust-builder

WORKDIR /app/aa56-rust
COPY aa56-rust/Cargo.toml ./Cargo.toml
COPY aa56-rust/catalogs ./catalogs
COPY aa56-rust/src ./src

RUN cargo build --release

# ==============================================================================
# Estágio 2: Runtime Leve Node.js + Binário Nativo Rust para Google Cloud Run
# ==============================================================================
FROM node:22-slim

WORKDIR /app

# Copia os catálogos astronômicos e o binário nativo Linux compilado no Estágio 1
COPY --from=rust-builder /app/aa56-rust/catalogs ./aa56-rust/catalogs
COPY --from=rust-builder /app/aa56-rust/target/release/aa56-cli ./aa56-rust/target/release/aa56-cli
RUN chmod +x ./aa56-rust/target/release/aa56-cli

# Copia o servidor e as páginas web da aplicação Node.js
COPY node-test-app ./node-test-app
COPY package.json ./package.json

# Google Cloud Run injeta automaticamente a variável PORT (padrão 8080)
ENV NODE_ENV=production
ENV PORT=8080
EXPOSE 8080

CMD ["node", "node-test-app/server.js"]

