FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_hail_gcloud:0.2.134.cpg1

# DeepVariant pa pipeline version.
ENV VERSION=0.1.2

# Add in the additional requirements that are most likely to change.
COPY LICENSE pyproject.toml README.md ./
COPY src src/
# COPY third_party third_party/

RUN apt-get update && apt-get install -y git

RUN pip install . \
    && pip install typing-extensions==4.12.0